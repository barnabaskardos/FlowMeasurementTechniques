from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "outputs" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "figure.dpi": 130,
        "savefig.dpi": 300,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "lines.linewidth": 1.6,
    }
)


COLORS = {
    0: "#1f77b4",
    5: "#2ca02c",
    15: "#d62728",
}


def read_table(name: str) -> pd.DataFrame:
    path = TABLES / name
    if not path.exists():
        raise FileNotFoundError(f"Run MATLAB processing first; missing {path}")
    return pd.read_csv(path)


def bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin({"true", "1"})


def save(fig: plt.Figure, name: str) -> None:
    if not fig.get_constrained_layout():
        fig.tight_layout()
    fig.savefig(FIGURES / f"{name}.png", bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_hwa_calibration() -> None:
    points = read_table("hwa_calibration_points.csv")
    curve = read_table("hwa_calibration_curve.csv")

    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.plot(curve["signal_v"], curve["velocity_pchip_ms"], label="PCHIP calibration")
    ax.plot(curve["signal_v"], curve["velocity_poly_ms"], "--", label="4th order polynomial")
    ax.scatter(points["signal_v"], points["velocity_ms"], s=22, color="black", zorder=3, label="Measured points")
    ax.set_xlabel("CTA output signal [V]")
    ax.set_ylabel("Velocity [m/s]")
    ax.set_title("Hot-wire calibration")
    ax.legend()
    save(fig, "hwa_calibration_curve")


def plot_hwa_profiles() -> None:
    profiles = read_table("hwa_profiles.csv")

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 4.2), sharey=True)
    for aoa, group in profiles.groupby("aoa_deg"):
        group = group.sort_values("y_relative_mm")
        color = COLORS.get(int(aoa), None)
        label = f"AoA {int(aoa)} deg"
        axes[0].plot(group["mean_velocity_ms"], group["y_relative_mm"], "-o", ms=3.5, color=color, label=label)
        axes[1].plot(group["rms_velocity_ms"], group["y_relative_mm"], "-o", ms=3.5, color=color, label=label)

    axes[0].set_xlabel("Mean velocity [m/s]")
    axes[0].set_ylabel("y relative to trailing edge [mm]")
    axes[0].set_title("HWA mean profiles")
    axes[1].set_xlabel("RMS fluctuation [m/s]")
    axes[1].set_title("HWA fluctuation profiles")
    axes[0].axhline(0, color="0.2", lw=0.8)
    axes[1].axhline(0, color="0.2", lw=0.8)
    axes[0].legend()
    save(fig, "hwa_velocity_profiles")


def plot_hwa_autocorrelation() -> None:
    autocorr = read_table("hwa_autocorrelation.csv")
    convergence = read_table("hwa_sampling_convergence.csv").iloc[0]
    visible = autocorr[autocorr["lag_s"] <= min(0.25, autocorr["lag_s"].max())]

    fig, ax = plt.subplots(figsize=(5.8, 3.6))
    ax.plot(visible["lag_s"], visible["rho"], color="#1f77b4")
    ax.axhline(0, color="black", lw=0.8)
    ax.axvline(convergence["first_zero_crossing_s"], color="#d62728", ls="--", label="First zero crossing")
    ax.set_xlabel("Lag [s]")
    ax.set_ylabel("Autocorrelation coefficient")
    ax.set_title("HWA autocorrelation used for sampling-time estimate")
    ax.legend()
    save(fig, "hwa_autocorrelation")


def plot_hwa_spectra() -> None:
    spectra = read_table("hwa_spectra.csv")
    peaks = read_table("hwa_spectral_peaks.csv")
    spectra = spectra[(spectra["frequency_hz"] >= 1) & (spectra["frequency_hz"] <= 5000)]

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    for (region, aoa), group in spectra.groupby(["region", "aoa_deg"]):
        color = COLORS.get(int(aoa), None)
        label = f"{region}, AoA {int(aoa)} deg"
        ax.loglog(group["frequency_hz"], group["phi_uu"], color=color, label=label)

    for _, peak in peaks.dropna(subset=["frequency_hz"]).iterrows():
        if 1 <= peak["frequency_hz"] <= 5000:
            ax.axvline(peak["frequency_hz"], color=COLORS.get(int(peak["aoa_deg"]), "0.5"), alpha=0.12)

    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Power spectral density [arb. units]")
    ax.set_title("HWA energy spectra")
    ax.legend()
    save(fig, "hwa_energy_spectra")


def field_grid(df: pd.DataFrame, value: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if "valid" in df.columns:
        valid = bool_series(df["valid"])
    else:
        valid = pd.Series(True, index=df.index)
    working = df.copy()
    working.loc[~valid, value] = np.nan
    if "v_normal_ms" in working.columns:
        working.loc[~valid, "v_normal_ms"] = np.nan

    value_grid = working.pivot(index="y_mm", columns="x_mm", values=value).sort_index().sort_index(axis=1)
    x = value_grid.columns.to_numpy(dtype=float)
    y = value_grid.index.to_numpy(dtype=float)
    X, Y = np.meshgrid(x, y)
    Z = value_grid.to_numpy(dtype=float)

    if "v_normal_ms" in working.columns and value != "v_normal_ms":
        v_grid = working.pivot(index="y_mm", columns="x_mm", values="v_normal_ms").sort_index().sort_index(axis=1)
        V = v_grid.to_numpy(dtype=float)
    else:
        V = np.zeros_like(Z)
    return X, Y, Z, V


def field_scale(df: pd.DataFrame, value: str = "u_streamwise_ms") -> tuple[float, float]:
    if "valid" in df.columns:
        valid = bool_series(df["valid"])
        values = df.loc[valid, value].to_numpy(dtype=float)
    else:
        values = df[value].to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return -1.0, 1.0
    vmin, vmax = np.nanpercentile(values, [2, 98])
    if np.isclose(vmin, vmax):
        pad = max(abs(vmin) * 0.1, 1.0)
        return vmin - pad, vmax + pad
    return float(vmin), float(vmax)


def draw_field_panel(ax: plt.Axes, df: pd.DataFrame, title: str, vmin: float, vmax: float):
    X, Y, U, V = field_grid(df, "u_streamwise_ms")
    levels = np.linspace(vmin, vmax, 25)
    contour = ax.contourf(X, Y, U, levels=levels, cmap="viridis", extend="both")
    skip_y = max(1, U.shape[0] // 16)
    skip_x = max(1, U.shape[1] // 18)
    ax.quiver(
        X[::skip_y, ::skip_x],
        Y[::skip_y, ::skip_x],
        U[::skip_y, ::skip_x],
        V[::skip_y, ::skip_x],
        color="white",
        scale=260,
        width=0.0024,
        alpha=0.85,
    )
    ax.set_title(title)
    ax.set_xlabel("x [mm]")
    ax.set_xticks([-100, -50, 0, 50])
    ax.set_aspect("equal", adjustable="box")
    return contour


def plot_piv_fields(table_name: str, figure_name: str, title_prefix: str) -> None:
    fields = read_table(table_name)
    fields = fields[fields["processing"] == "Overlap50MP3"]
    aoa_values = sorted(fields["aoa_deg"].unique())
    fig, axes = plt.subplots(
        1,
        len(aoa_values),
        figsize=(12.2, 3.7),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    if len(aoa_values) == 1:
        axes = [axes]

    for ax, aoa in zip(axes, aoa_values):
        subset = fields[fields["aoa_deg"] == aoa]
        X, Y, U, V = field_grid(subset, "u_streamwise_ms")
        contour = ax.contourf(X, Y, U, levels=24, cmap="viridis")
        skip_y = max(1, U.shape[0] // 16)
        skip_x = max(1, U.shape[1] // 18)
        ax.quiver(X[::skip_y, ::skip_x], Y[::skip_y, ::skip_x], U[::skip_y, ::skip_x], V[::skip_y, ::skip_x],
                  color="white", scale=260, width=0.0024, alpha=0.85)
        ax.axhline(55, color="white", lw=0.8, alpha=0.65)
        ax.set_title(f"{title_prefix}, AoA {int(aoa)} deg")
        ax.set_xlabel("x [mm]")
        ax.set_xticks([-100, -50, 0, 50])
        ax.set_aspect("equal", adjustable="box")
    axes[0].set_ylabel("y [mm]")
    cbar = fig.colorbar(contour, ax=axes, shrink=0.86)
    cbar.set_label("Streamwise velocity [m/s]")
    save(fig, figure_name)


def plot_piv_profile_comparison() -> None:
    comparison = read_table("technique_profile_comparison.csv")
    comparison = comparison[bool_series(comparison["valid"])]
    aoa_values = sorted(comparison["aoa_deg"].unique())

    fig, axes = plt.subplots(2, len(aoa_values), figsize=(12.0, 6.4), sharey=True)
    if len(aoa_values) == 1:
        axes = np.array([[axes[0]], [axes[1]]])

    markers = {"HWA": "o", "PIV": "s", "Pressure probe": "^"}
    for col, aoa in enumerate(aoa_values):
        subset = comparison[comparison["aoa_deg"] == aoa]
        for technique, group in subset.groupby("technique"):
            group = group.sort_values("y_relative_mm")
            marker = markers.get(technique, "o")
            axes[0, col].plot(group["mean_velocity_ms"], group["y_relative_mm"], marker=marker, ms=3.5, label=technique)
            if not group["rms_velocity_ms"].isna().all():
                axes[1, col].plot(group["rms_velocity_ms"], group["y_relative_mm"], marker=marker, ms=3.5, label=technique)
        axes[0, col].set_title(f"AoA {int(aoa)} deg")
        axes[1, col].set_xlabel("Velocity [m/s]")
        axes[0, col].axhline(0, color="0.2", lw=0.8)
        axes[1, col].axhline(0, color="0.2", lw=0.8)
    axes[0, 0].set_ylabel("y relative to TE [mm]")
    axes[1, 0].set_ylabel("y relative to TE [mm]")
    axes[0, 0].legend()
    axes[0, 0].text(0.02, 1.05, "Mean velocity", transform=axes[0, 0].transAxes, weight="bold")
    axes[1, 0].text(0.02, 1.05, "RMS fluctuation", transform=axes[1, 0].transAxes, weight="bold")
    save(fig, "technique_profile_comparison")


def plot_piv_processing_summary() -> None:
    summary = read_table("piv_processing_summary.csv")
    interesting = summary[
        summary["processing"].str.contains("Overlap50SinglePass|Overlap50MP3|10Img", regex=True)
    ].copy()
    if interesting.empty:
        return
    interesting["label"] = interesting["case_id"].astype(str) + "\n" + interesting["processing"].astype(str)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    x = np.arange(len(interesting))
    axes[0].bar(x, interesting["valid_fraction"], color="#4c78a8")
    axes[0].set_ylabel("Valid-vector fraction")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(interesting["label"], rotation=70, ha="right")
    axes[0].set_ylim(0, 1.05)

    axes[1].bar(x, interesting["spatial_mean_u_rms_ms"], color="#f58518")
    axes[1].set_ylabel("Spatial mean u RMS [m/s]")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(interesting["label"], rotation=70, ha="right")
    axes[1].set_title("Processing and ensemble-size sensitivity")
    save(fig, "piv_processing_summary")


def plot_davis_window_size_screenshots() -> None:
    screenshot_dir = ROOT / "data" / "PIV processed" / "Screenshots"
    files = [
        ("16 x 16 px", screenshot_dir / "16x16Test.png"),
        ("32 x 32 px", screenshot_dir / "32x32Test.png"),
        ("64 x 64 px", screenshot_dir / "64x64Test.png"),
    ]
    files = [(label, path) for label, path in files if path.exists()]
    if not files:
        return

    fig, axes = plt.subplots(1, len(files), figsize=(12.2, 4.2), constrained_layout=True)
    if len(files) == 1:
        axes = [axes]
    for ax, (label, path) in zip(axes, files):
        image = plt.imread(path)
        ax.imshow(image)
        ax.set_title(f"DaVis window-size test\n{label}")
        ax.axis("off")
    save(fig, "piv_davis_window_size_screenshots")


def plot_self_piv_window_fields() -> None:
    fields = read_table("piv_self_fields.csv")
    windows = [16, 32, 64]
    selected = fields[fields["window_size_px"].isin(windows)]
    if selected.empty:
        return
    vmin, vmax = field_scale(selected)

    fig, axes = plt.subplots(1, len(windows), figsize=(12.2, 4.0), sharex=False, sharey=False, constrained_layout=True)
    for ax, window in zip(axes, windows):
        subset = selected[selected["window_size_px"] == window]
        contour = draw_field_panel(ax, subset, f"Self-made PIV\n{window} x {window} px", vmin, vmax)
    axes[0].set_ylabel("y [mm]")
    cbar = fig.colorbar(contour, ax=axes, shrink=0.86)
    cbar.set_label("Streamwise velocity [m/s]")
    save(fig, "piv_self_window_fields")


def plot_parameter_field_group(group: str, labels: list[str], figure_name: str, title: str) -> None:
    fields = read_table("piv_parameter_fields.csv")
    selected = fields[fields["parameter_group"] == group]
    selected = selected[selected["parameter_label"].isin(labels)]
    if selected.empty:
        return

    present_labels = [label for label in labels if label in set(selected["parameter_label"])]
    vmin, vmax = field_scale(selected)
    fig, axes = plt.subplots(1, len(present_labels), figsize=(12.2, 4.0), sharex=True, sharey=True, constrained_layout=True)
    if len(present_labels) == 1:
        axes = [axes]

    for ax, label in zip(axes, present_labels):
        subset = selected[selected["parameter_label"] == label]
        contour = draw_field_panel(ax, subset, label, vmin, vmax)
    axes[0].set_ylabel("y [mm]")
    fig.suptitle(title)
    cbar = fig.colorbar(contour, ax=axes, shrink=0.86)
    cbar.set_label("Streamwise velocity [m/s]")
    save(fig, figure_name)


def plot_piv_parameter_studies() -> None:
    plot_parameter_field_group(
        "overlap_and_pass",
        ["0% overlap, single pass", "50% overlap, single pass", "50% overlap, 3-pass"],
        "piv_overlap_multipass_comparison",
        "Effect of overlap and single-pass vs multipass processing",
    )
    plot_parameter_field_group(
        "ensemble_size",
        ["10 sequential samples", "10 separated samples", "100 samples"],
        "piv_ensemble_size_comparison",
        "Effect of ensemble size and sample spacing on the mean field",
    )
    plot_parameter_field_group(
        "delta_t",
        ["Original Delta t", "Short Delta t"],
        "piv_delta_t_comparison",
        "Effect of pulse separation on the processed field",
    )


def plot_piv_frame_brightness() -> None:
    image_path = ROOT / "data" / "PIV" / "FMT Results" / "aoa15_final" / "B00001.tif"
    if not image_path.exists():
        return

    image = plt.imread(image_path)
    if image.ndim == 3:
        image = image.mean(axis=2)
    half = image.shape[0] // 2
    upper = image[:half, :]
    lower = image[half : 2 * half, :]
    combined = np.concatenate([upper.ravel(), lower.ravel()])
    vmin, vmax = np.nanpercentile(combined, [1, 99.7])

    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.2), constrained_layout=True)
    axes[0, 0].imshow(upper, cmap="gray", vmin=vmin, vmax=vmax)
    axes[0, 0].set_title(f"Upper frame, mean intensity {upper.mean():.0f}")
    axes[0, 1].imshow(lower, cmap="gray", vmin=vmin, vmax=vmax)
    axes[0, 1].set_title(f"Lower frame, mean intensity {lower.mean():.0f}")
    for ax in axes[0, :]:
        ax.axis("off")

    bins = np.linspace(vmin, vmax, 120)
    axes[1, 0].hist(upper.ravel(), bins=bins, density=True, alpha=0.75, label="Upper frame")
    axes[1, 0].hist(lower.ravel(), bins=bins, density=True, alpha=0.65, label="Lower frame")
    axes[1, 0].set_xlabel("Pixel intensity")
    axes[1, 0].set_ylabel("Probability density")
    axes[1, 0].set_title("Intensity distributions")
    axes[1, 0].legend()

    difference = lower.astype(float) - upper.astype(float)
    axes[1, 1].imshow(difference, cmap="coolwarm", vmin=-np.nanpercentile(abs(difference), 99), vmax=np.nanpercentile(abs(difference), 99))
    axes[1, 1].set_title("Lower - upper intensity")
    axes[1, 1].axis("off")
    save(fig, "piv_frame_brightness")


def plot_self_piv() -> None:
    sensitivity = read_table("piv_self_window_sensitivity.csv")
    comparison = read_table("piv_self_vs_davis_points.csv")

    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    ax.plot(sensitivity["window_size_px"], sensitivity["mean_speed_ms"], "-o", label="Mean speed")
    ax.plot(sensitivity["window_size_px"], sensitivity["spatial_std_u_streamwise_ms"], "-s", label="Spatial std(u)")
    ax.set_xlabel("Interrogation window size [px]")
    ax.set_ylabel("Velocity scale [m/s]")
    ax.set_title("Self-made PIV window-size sensitivity")
    ax.legend()
    save(fig, "piv_self_window_sensitivity")

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.7), sharex=True, sharey=True, constrained_layout=True)
    for ax, value, title in zip(
        axes,
        ["davis_u_streamwise_ms", "self_u_streamwise_ms", "difference_u_ms"],
        ["DaVis u", "Self-made u", "Self - DaVis"],
    ):
        X, Y, Z, _ = field_grid(comparison.rename(columns={value: "u_streamwise_ms"}), "u_streamwise_ms")
        cmap = "coolwarm" if "difference" in value else "viridis"
        levels = 25
        contour = ax.contourf(X, Y, Z, levels=levels, cmap=cmap)
        ax.set_title(title)
        ax.set_xlabel("x [mm]")
        ax.set_aspect("equal", adjustable="box")
        fig.colorbar(contour, ax=ax, shrink=0.82)
    axes[0].set_ylabel("y [mm]")
    save(fig, "piv_self_vs_davis")


def plot_frequency_resolution() -> None:
    piv = read_table("piv_frequency_resolution.csv")
    hwa = read_table("hwa_frequency_resolution.csv")

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6))
    axes[0].bar(piv["ensemble"], piv["nyquist_hz"], color="#4c78a8")
    axes[0].set_ylabel("Nyquist frequency [Hz]")
    axes[0].set_title("PIV resolvable frequency range")
    axes[0].tick_params(axis="x", rotation=20)

    labels = hwa["region"].astype(str) + "\nAoA " + hwa["aoa_deg"].astype(int).astype(str)
    axes[1].bar(labels, hwa["frequency_resolution_hz"], color="#f58518")
    axes[1].set_ylabel("Frequency resolution [Hz]")
    axes[1].set_title("HWA Welch resolution")
    axes[1].tick_params(axis="x", rotation=20)
    save(fig, "frequency_resolution")


def main() -> None:
    plot_hwa_calibration()
    plot_hwa_profiles()
    plot_hwa_autocorrelation()
    plot_hwa_spectra()
    plot_piv_fields("piv_mean_fields.csv", "piv_mean_fields", "Mean PIV field")
    plot_piv_fields("piv_instantaneous_fields.csv", "piv_instantaneous_fields", "Instantaneous PIV field")
    plot_piv_profile_comparison()
    plot_piv_processing_summary()
    plot_davis_window_size_screenshots()
    plot_self_piv_window_fields()
    plot_piv_parameter_studies()
    plot_piv_frame_brightness()
    plot_self_piv()
    plot_frequency_resolution()
    print(f"Figures written to {FIGURES}")


if __name__ == "__main__":
    main()
