# Report Notes For PIV Parameter Discussion

Use these points as working notes for the Results and Discussion section.

## PIV Processing Parameters

Recommended figures:

- `outputs/figures/piv_davis_window_size_screenshots.png`
- `outputs/figures/piv_self_window_fields.png`
- `outputs/figures/piv_self_window_sensitivity.png`
- `outputs/figures/piv_overlap_multipass_comparison.png`
- `outputs/figures/piv_ensemble_size_comparison.png`
- `outputs/figures/piv_delta_t_comparison.png`
- `outputs/figures/piv_frame_brightness.png`

Why 32 x 32 px was selected:

- The 16 x 16 px result has the highest vector density, but the field is visibly
  noisy and contains many small-scale velocity fluctuations that are not
  physically coherent. In the self-made PIV result, the spatial standard
  deviation of streamwise velocity is highest for this window size.
- The 64 x 64 px result is smoother, but it reduces spatial resolution and
  smears the wake/shear-layer gradients near the airfoil.
- The 32 x 32 px result is the compromise: it keeps enough spatial resolution to
  resolve the wake and separated region while strongly reducing the noisy
  small-scale behavior seen with 16 x 16 px.
- The self-made PIV summary table supports this: `piv_self_window_sensitivity.csv`
  shows high spatial scatter for 16 px, while 32 px and 64 px are similar in
  stability. The 32 px case is preferred because it preserves twice the spatial
  resolution of 64 px in each direction.

Overlap and single-pass/multipass:

- 50% overlap increases vector density compared with 0% overlap, making the wake
  boundary and shear layer easier to identify.
- The 3-pass result is close to the 50% overlap single-pass result in the
  current exported data, but it gives a slightly smoother and more coherent
  deformation-aware field. This is the standard professional-processing choice.

Ensemble size:

- The 10-image averages recover the main wake pattern, but the 100-image average
  is the better report result because random fluctuations are reduced and the
  mean wake topology is more stable.
- The 10 separated samples and 10 sequential samples can be compared to discuss
  whether sample spacing changes the mean. In these data the large-scale
  structure is similar, so the main benefit of 100 samples is convergence rather
  than a changed mean topology.

Pulse separation:

- The short-Delta-t case reduces particle displacement, which can help avoid
  excessive displacements in high-speed regions, but it also makes the measured
  displacement smaller relative to pixel and correlation uncertainty.
- The original Delta-t case is therefore more suitable for the main velocity
  field when the one-quarter-window displacement rule is satisfied.

Frame brightness:

- `piv_frame_brightness.png` shows that the two exported frames do not have the
  same background/intensity distribution. Discuss this as a double-frame camera
  exposure effect and as a reason why background subtraction and normalized
  correlation are useful.

Self-made PIV vs DaVis:

- `piv_self_vs_davis.png` and `piv_self_vs_davis_summary.csv` show the expected
  accuracy gap between a simple direct cross-correlation implementation and
  professional DaVis processing.
- The self-made code uses fixed windows and no overlap, no subpixel peak fit, no
  multipass window shifting/deformation, and no vector validation/replacement.
  These omissions explain most of the difference in accuracy.
