clear;
clc;
close all;

cfg = defaultConfig();
ensureOutputFolders(cfg);

fprintf('Flow Measurement Techniques processing pipeline\n');
fprintf('Repository root: %s\n', cfg.repoRoot);
fprintf('Output folder:   %s\n\n', cfg.outputDir);

[calibration, calRaw, calFit, calPoly] = loadHwaCalibration(cfg);
writetable(calRaw, fullfile(cfg.tablesDir, 'hwa_calibration_points.csv'));
writetable(calFit, fullfile(cfg.tablesDir, 'hwa_calibration_curve.csv'));
writetable(calPoly, fullfile(cfg.tablesDir, 'hwa_calibration_polynomial.csv'));

[hwaProfiles, convergence, autocorrelation, spectra, spectralPeaks, hwaFrequency] = ...
    processHwaData(cfg, calibration);
writetable(hwaProfiles, fullfile(cfg.tablesDir, 'hwa_profiles.csv'));
writetable(convergence, fullfile(cfg.tablesDir, 'hwa_sampling_convergence.csv'));
writetable(autocorrelation, fullfile(cfg.tablesDir, 'hwa_autocorrelation.csv'));
writetable(spectra, fullfile(cfg.tablesDir, 'hwa_spectra.csv'));
writetable(spectralPeaks, fullfile(cfg.tablesDir, 'hwa_spectral_peaks.csv'));
writetable(hwaFrequency, fullfile(cfg.tablesDir, 'hwa_frequency_resolution.csv'));

[pivMeanFields, pivRmsFields, pivInstantFields, pivProfiles, ...
    pivProcessingSummary, pivFrequency, selfFields, selfWindowSummary, ...
    selfComparePoints, selfCompareSummary, pivParameterFields, ...
    pivParameterSummary] = processPivData(cfg);
writetable(pivMeanFields, fullfile(cfg.tablesDir, 'piv_mean_fields.csv'));
writetable(pivRmsFields, fullfile(cfg.tablesDir, 'piv_rms_fields.csv'));
writetable(pivInstantFields, fullfile(cfg.tablesDir, 'piv_instantaneous_fields.csv'));
writetable(pivProfiles, fullfile(cfg.tablesDir, 'piv_profile_xc12.csv'));
writetable(pivProcessingSummary, fullfile(cfg.tablesDir, 'piv_processing_summary.csv'));
writetable(pivFrequency, fullfile(cfg.tablesDir, 'piv_frequency_resolution.csv'));
writetable(selfFields, fullfile(cfg.tablesDir, 'piv_self_fields.csv'));
writetable(selfWindowSummary, fullfile(cfg.tablesDir, 'piv_self_window_sensitivity.csv'));
writetable(selfComparePoints, fullfile(cfg.tablesDir, 'piv_self_vs_davis_points.csv'));
writetable(selfCompareSummary, fullfile(cfg.tablesDir, 'piv_self_vs_davis_summary.csv'));
writetable(pivParameterFields, fullfile(cfg.tablesDir, 'piv_parameter_fields.csv'));
writetable(pivParameterSummary, fullfile(cfg.tablesDir, 'piv_parameter_summary.csv'));

comparisonProfile = buildTechniqueComparisonProfile(cfg, hwaProfiles, pivProfiles);
writetable(comparisonProfile, fullfile(cfg.tablesDir, 'technique_profile_comparison.csv'));

missingInputs = detectMissingInputs(cfg);
writetable(missingInputs, fullfile(cfg.tablesDir, 'missing_required_inputs.csv'));

fprintf('\nDone. MATLAB processing tables were written to:\n%s\n', cfg.tablesDir);
fprintf('Next step: python analysis/plot_results.py\n');

function cfg = defaultConfig()
    scriptDir = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(scriptDir);
    if ~isfolder(fullfile(repoRoot, 'data'))
        repoRoot = scriptDir;
    end

    cfg.repoRoot = repoRoot;
    cfg.dataDir = fullfile(repoRoot, 'data');
    cfg.outputDir = fullfile(repoRoot, 'outputs');
    cfg.tablesDir = fullfile(cfg.outputDir, 'tables');
    cfg.figuresDir = fullfile(cfg.outputDir, 'figures');

    cfg.hwaDir = fullfile(cfg.dataDir, 'HWA', 'Group15');
    cfg.pivRawDir = fullfile(cfg.dataDir, 'PIV', 'FMT Results');
    cfg.pivProcessedDir = fullfile(cfg.dataDir, 'PIV processed');

    cfg.trailingEdgePosition_mm = 55;
    cfg.pivProfileX_mm = 20;
    cfg.pivStreamwiseSign = -1;
    cfg.pivSampleRate_Hz = 8.3;
    cfg.hwaTargetUncertainty = 0.01;
    cfg.hwaCoverageFactor = 3;
    cfg.hwaAutocorrMaxLag_s = 1;
    cfg.spectraWindowSamples = 16384;
    cfg.spectraOverlapFraction = 0.5;
    cfg.spectraPlotMaxFrequency_Hz = 5000;
    cfg.minimumPeakFrequency_Hz = 10;
    cfg.minimumPeakSeparation_Hz = 20;
    cfg.nDominantPeaks = 3;

    cfg.selfPiv.windowSizes_px = [16, 32, 64];
    cfg.selfPiv.deltaT_s = 74e-6;
    cfg.selfPiv.imagePath = fullfile(cfg.pivRawDir, 'aoa15_final', 'B00001.tif');
    cfg.selfPiv.calibrationImagePath = fullfile(cfg.pivRawDir, 'cal_final', 'B00001.tif');
    cfg.selfPiv.referenceDavisPath = fullfile(cfg.pivProcessedDir, ...
        'Processed15', 'Overlap0SinglePass', 'B00001.dat');

    cfg.pivCases = struct( ...
        'aoa_deg', {0, 5, 15}, ...
        'case_id', {'aoa0', 'aoa5', 'aoa15'}, ...
        'processedFolder', {'ProcessedAoa0', 'Processed5', 'Processed15'}, ...
        'rawImagePath', { ...
            fullfile(cfg.pivRawDir, 'AoA0_final', 'B00001.tif'), ...
            fullfile(cfg.pivRawDir, 'aoa5_final', 'B00001.tif'), ...
            fullfile(cfg.pivRawDir, 'aoa15_final', 'B00001.tif')}, ...
        'nominalSamples', {20, 20, 100});

    cfg.pivShortDtCase = struct( ...
        'aoa_deg', 15, ...
        'case_id', 'aoa15_short_dt', ...
        'processedFolder', 'Processed15dt', ...
        'rawImagePath', fullfile(cfg.pivRawDir, 'aoa15deltat_final', 'B00001.tif'), ...
        'nominalSamples', 20);
end

function ensureOutputFolders(cfg)
    folders = {cfg.outputDir, cfg.tablesDir, cfg.figuresDir};
    for iFolder = 1:numel(folders)
        if ~isfolder(folders{iFolder})
            mkdir(folders{iFolder});
        end
    end
end

function [calibration, rawTable, fitTable, polyTable] = loadHwaCalibration(cfg)
    calibrationPath = fullfile(cfg.hwaDir, 'calibration.txt');
    data = readmatrix(calibrationPath);
    data = data(:, 1:2);
    data = data(all(isfinite(data), 2), :);

    velocity = data(:, 1);
    signal = data(:, 2);
    [signalUnique, ~, groupIdx] = unique(signal, 'sorted');
    velocityUnique = accumarray(groupIdx, velocity, [], @mean);

    polyOrder = min(4, numel(signalUnique) - 1);
    polyCoefficients = polyfit(signalUnique, velocityUnique, polyOrder);

    signalGrid = linspace(min(signalUnique), max(signalUnique), 400).';
    fitTable = table( ...
        signalGrid, ...
        interp1(signalUnique, velocityUnique, signalGrid, 'pchip'), ...
        polyval(polyCoefficients, signalGrid), ...
        'VariableNames', {'signal_v', 'velocity_pchip_ms', 'velocity_poly_ms'});

    rawTable = table(velocity, signal, ...
        'VariableNames', {'velocity_ms', 'signal_v'});

    powers = polyOrder:-1:0;
    polyTable = table(powers.', polyCoefficients.', ...
        'VariableNames', {'power', 'coefficient'});

    calibration.signalUnique = signalUnique;
    calibration.velocityUnique = velocityUnique;
    calibration.polyCoefficients = polyCoefficients;
end

function velocity = convertHwaSignal(calibration, signal)
    velocity = interp1(calibration.signalUnique, calibration.velocityUnique, ...
        signal, 'pchip', 'extrap');
end

function [profiles, convergence, autocorrelation, spectra, spectralPeaks, frequencyResolution] = ...
        processHwaData(cfg, calibration)
    files = dir(fullfile(cfg.hwaDir, 'aoa*_10.txt'));
    records = struct('aoa_deg', {}, 'position_mm', {}, 'y_relative_mm', {}, ...
        'mean_velocity_ms', {}, 'rms_velocity_ms', {}, 'sample_rate_hz', {}, ...
        'n_samples', {}, 'file_name', {});

    for iFile = 1:numel(files)
        fileName = files(iFile).name;
        tokens = regexp(fileName, '^aoa(-?\d+)_([0-9]+)_10\.txt$', 'tokens', 'once');
        if isempty(tokens)
            continue;
        end

        aoa = str2double(tokens{1});
        position = str2double(tokens{2});
        filePath = fullfile(files(iFile).folder, fileName);
        [time, velocity] = readHwaVelocity(filePath, calibration);

        records(end + 1).aoa_deg = aoa; %#ok<AGROW>
        records(end).position_mm = position;
        records(end).y_relative_mm = position - cfg.trailingEdgePosition_mm;
        records(end).mean_velocity_ms = mean(velocity, 'omitnan');
        records(end).rms_velocity_ms = rmsAboutMean(velocity);
        records(end).sample_rate_hz = 1 / mean(diff(time));
        records(end).n_samples = numel(velocity);
        records(end).file_name = string(fileName);
    end

    if isempty(records)
        error('No HWA measurement records were found in %s.', cfg.hwaDir);
    end

    profiles = sortrows(struct2table(records), {'aoa_deg', 'position_mm'});

    [convergence, autocorrelation] = computeSamplingConvergence(cfg, calibration);
    [spectra, spectralPeaks, frequencyResolution] = computeHwaSpectra(cfg, calibration, profiles);
end

function [time, velocity] = readHwaVelocity(filePath, calibration)
    rawData = readmatrix(filePath);
    time = rawData(:, 1);
    signal = rawData(:, 2);
    valid = isfinite(time) & isfinite(signal);
    time = time(valid);
    signal = signal(valid);
    velocity = convertHwaSignal(calibration, signal);
end

function value = rmsAboutMean(signal)
    meanValue = mean(signal, 'omitnan');
    value = sqrt(mean((signal - meanValue).^2, 'omitnan'));
end

function [convergence, autocorrelation] = computeSamplingConvergence(cfg, calibration)
    filePath = fullfile(cfg.hwaDir, 'correlationtest.txt');
    [time, velocity] = readHwaVelocity(filePath, calibration);
    sampleRate = 1 / mean(diff(time));
    fluctuation = velocity - mean(velocity, 'omitnan');
    maxLagSamples = min(numel(fluctuation) - 1, round(cfg.hwaAutocorrMaxLag_s * sampleRate));

    [rhoFull, lags] = xcorr(fluctuation, maxLagSamples, 'coeff');
    positive = lags >= 0;
    lag_s = lags(positive);
    lag_s = lag_s(:) / sampleRate;
    rho = rhoFull(positive);
    rho = rho(:);

    firstZeroIndex = find(rho <= 0 & lag_s > 0, 1, 'first');
    if isempty(firstZeroIndex)
        firstZeroIndex = numel(rho);
    end

    integrationIndex = 1:firstZeroIndex;
    integralTime = trapz(lag_s(integrationIndex), rho(integrationIndex));
    meanVelocity = mean(velocity, 'omitnan');
    rmsVelocity = rmsAboutMean(velocity);
    requiredSamples = (cfg.hwaCoverageFactor * rmsVelocity / ...
        (cfg.hwaTargetUncertainty * abs(meanVelocity)))^2;
    requiredSamplingTime = 2 * integralTime * requiredSamples;

    convergence = table( ...
        meanVelocity, rmsVelocity, sampleRate, lag_s(firstZeroIndex), ...
        integralTime, cfg.hwaCoverageFactor, cfg.hwaTargetUncertainty, ...
        requiredSamples, requiredSamplingTime, ...
        'VariableNames', {'mean_velocity_ms', 'rms_velocity_ms', ...
        'sample_rate_hz', 'first_zero_crossing_s', 'integral_time_s', ...
        'coverage_factor', 'target_uncertainty_fraction', ...
        'uncorrelated_samples_required', 'sampling_time_required_s'});

    autocorrelation = table(lag_s, rho, ...
        'VariableNames', {'lag_s', 'rho'});
end

function [spectra, spectralPeaks, frequencyResolution] = computeHwaSpectra(cfg, calibration, profiles)
    aoa15 = profiles(profiles.aoa_deg == 15, :);
    [~, shearIndex] = max(aoa15.rms_velocity_ms);
    shearLayerPosition = aoa15.position_mm(shearIndex);

    requests = table( ...
        ["Trailing edge"; "Trailing edge"; "Shear layer"], ...
        [cfg.trailingEdgePosition_mm; cfg.trailingEdgePosition_mm; shearLayerPosition], ...
        [0; 5; 15], ...
        'VariableNames', {'region', 'position_mm', 'aoa_deg'});

    spectra = table();
    spectralPeaks = table();
    frequencyResolution = table();

    for iRequest = 1:height(requests)
        aoa = requests.aoa_deg(iRequest);
        position = requests.position_mm(iRequest);
        region = requests.region(iRequest);
        fileName = sprintf('aoa%d_%d_10.txt', aoa, position);
        filePath = fullfile(cfg.hwaDir, fileName);

        [time, velocity] = readHwaVelocity(filePath, calibration);
        [frequency_hz, phi_uu, df_hz, sampleRate] = computeWelchSpectrum( ...
            time, velocity, cfg.spectraWindowSamples, cfg.spectraOverlapFraction);

        thisSpectra = table( ...
            repmat(region, numel(frequency_hz), 1), ...
            repmat(position, numel(frequency_hz), 1), ...
            repmat(aoa, numel(frequency_hz), 1), ...
            frequency_hz(:), phi_uu(:), ...
            'VariableNames', {'region', 'position_mm', 'aoa_deg', ...
            'frequency_hz', 'phi_uu'});
        spectra = [spectra; thisSpectra]; %#ok<AGROW>

        peakFrequencies = dominantFrequencies(frequency_hz, phi_uu, ...
            cfg.minimumPeakFrequency_Hz, cfg.spectraPlotMaxFrequency_Hz, ...
            cfg.minimumPeakSeparation_Hz, cfg.nDominantPeaks);
        peakFrequencies(end + 1:cfg.nDominantPeaks) = NaN;
        peakNumbers = (1:cfg.nDominantPeaks).';
        peakTable = table( ...
            repmat(region, cfg.nDominantPeaks, 1), ...
            repmat(position, cfg.nDominantPeaks, 1), ...
            repmat(aoa, cfg.nDominantPeaks, 1), ...
            peakNumbers, peakFrequencies(:), ...
            'VariableNames', {'region', 'position_mm', 'aoa_deg', ...
            'peak_number', 'frequency_hz'});
        spectralPeaks = [spectralPeaks; peakTable]; %#ok<AGROW>

        resolutionTable = table(region, position, aoa, sampleRate, ...
            sampleRate / 2, df_hz, cfg.spectraWindowSamples, ...
            'VariableNames', {'region', 'position_mm', 'aoa_deg', ...
            'sample_rate_hz', 'nyquist_hz', 'frequency_resolution_hz', ...
            'window_samples'});
        frequencyResolution = [frequencyResolution; resolutionTable]; %#ok<AGROW>
    end
end

function [frequency_hz, phi_uu, df_hz, sampleRate] = computeWelchSpectrum( ...
        time, velocity, windowSamples, overlapFraction)
    time = time(:);
    velocity = velocity(:);
    velocity = velocity - mean(velocity, 'omitnan');
    sampleRate = 1 / mean(diff(time));
    windowSamples = min(windowSamples, numel(velocity));
    overlapSamples = floor(overlapFraction * windowSamples);
    window = hann(windowSamples, 'periodic');
    [phi_uu, frequency_hz] = pwelch(velocity, window, overlapSamples, ...
        windowSamples, sampleRate);
    df_hz = sampleRate / windowSamples;
end

function peakFrequencies_hz = dominantFrequencies(frequency_hz, phi_uu, ...
        minFrequency_hz, maxFrequency_hz, minSeparation_hz, nPeaks)
    valid = frequency_hz >= minFrequency_hz & frequency_hz <= maxFrequency_hz;
    candidate = false(size(valid));
    candidate(2:end - 1) = valid(2:end - 1) & ...
        phi_uu(2:end - 1) > phi_uu(1:end - 2) & ...
        phi_uu(2:end - 1) >= phi_uu(3:end);

    peakFrequencies = frequency_hz(candidate);
    peakEnergies = phi_uu(candidate);
    [~, order] = sort(peakEnergies, 'descend');
    peakFrequencies = peakFrequencies(order);

    peakFrequencies_hz = [];
    for iPeak = 1:numel(peakFrequencies)
        if isempty(peakFrequencies_hz) || ...
                all(abs(peakFrequencies(iPeak) - peakFrequencies_hz) >= minSeparation_hz)
            peakFrequencies_hz(end + 1) = peakFrequencies(iPeak); %#ok<AGROW>
        end
        if numel(peakFrequencies_hz) == nPeaks
            break;
        end
    end
end

function [meanFields, rmsFields, instantFields, profiles, processingSummary, ...
        pivFrequency, selfFields, selfWindowSummary, selfComparePoints, ...
        selfCompareSummary, parameterFields, parameterSummary] = ...
        processPivData(cfg)
    meanFields = table();
    rmsFields = table();
    instantFields = table();
    profiles = table();
    processingSummary = table();

    for iCase = 1:numel(cfg.pivCases)
        thisCase = cfg.pivCases(iCase);
        caseDir = fullfile(cfg.pivProcessedDir, thisCase.processedFolder);

        [avgFile, stdevFile] = findAverageAndStdevFiles( ...
            fullfile(caseDir, 'Overlap50MP3AvgStDev'));
        avgField = readDavisField(avgFile, true);
        stdevField = readDavisField(stdevFile, true);

        meanFields = [meanFields; davisFieldToTable(avgField, thisCase, ...
            'Overlap50MP3', 'mean', cfg)]; %#ok<AGROW>
        rmsFields = [rmsFields; davisRmsFieldToTable(stdevField, thisCase, ...
            'Overlap50MP3', cfg)]; %#ok<AGROW>

        instantPath = fullfile(caseDir, 'Overlap50MP3', 'B00001.dat');
        if isfile(instantPath)
            instantField = readDavisField(instantPath, true);
            instantFields = [instantFields; davisFieldToTable(instantField, ...
                thisCase, 'Overlap50MP3', 'instantaneous', cfg)]; %#ok<AGROW>
        end

        profiles = [profiles; extractPivProfile(cfg, thisCase, ...
            avgField, stdevField, 'Overlap50MP3')]; %#ok<AGROW>

        processingSummary = [processingSummary; summarizePivProcessingFolder(cfg, thisCase)]; %#ok<AGROW>
    end

    processingSummary = [processingSummary; summarizePivProcessingFolder(cfg, cfg.pivShortDtCase)]; %#ok<AGROW>

    pivFrequency = table( ...
        ["20 snapshots"; "100 snapshots"], ...
        [20; 100], ...
        repmat(cfg.pivSampleRate_Hz, 2, 1), ...
        repmat(cfg.pivSampleRate_Hz / 2, 2, 1), ...
        cfg.pivSampleRate_Hz ./ [20; 100], ...
        [20; 100] ./ cfg.pivSampleRate_Hz, ...
        'VariableNames', {'ensemble', 'n_snapshots', 'sample_rate_hz', ...
        'nyquist_hz', 'frequency_resolution_hz', 'record_duration_s'});

    [selfFields, selfWindowSummary, selfComparePoints, selfCompareSummary] = processSelfPiv(cfg);
    [parameterFields, parameterSummary] = buildPivParameterStudy(cfg);
end

function [parameterFields, parameterSummary] = buildPivParameterStudy(cfg)
    aoa15Case = cfg.pivCases([cfg.pivCases.aoa_deg] == 15);
    specs = {
        "overlap_and_pass", "0% overlap, single pass", aoa15Case, "Overlap0SinglePassAvgStDev", 100, "original";
        "overlap_and_pass", "50% overlap, single pass", aoa15Case, "Overlap50SinglePassAvgStDev", 100, "original";
        "overlap_and_pass", "50% overlap, 3-pass", aoa15Case, "Overlap50MP3AvgStDev", 100, "original";
        "ensemble_size", "10 sequential samples", aoa15Case, "Overlap50MP3Img10Inc1AvgStDev", 10, "original";
        "ensemble_size", "10 separated samples", aoa15Case, "Overlap50MP3Img10Inc10AvgStDev", 10, "original";
        "ensemble_size", "100 samples", aoa15Case, "Overlap50MP3AvgStDev", 100, "original";
        "delta_t", "Original Delta t", aoa15Case, "Overlap50MP3AvgStDev", 100, "original";
        "delta_t", "Short Delta t", cfg.pivShortDtCase, "Overlap50MP3AvgStDev", 20, "short";
    };

    parameterFields = table();
    parameterSummary = table();

    for iSpec = 1:size(specs, 1)
        parameterGroup = specs{iSpec, 1};
        parameterLabel = specs{iSpec, 2};
        pivCase = specs{iSpec, 3};
        avgFolder = specs{iSpec, 4};
        nSnapshots = specs{iSpec, 5};
        deltaTLabel = specs{iSpec, 6};

        folderPath = fullfile(cfg.pivProcessedDir, pivCase.processedFolder, avgFolder);
        if ~isfolder(folderPath)
            warning('Skipping missing PIV parameter-study folder: %s', folderPath);
            continue;
        end

        [avgFile, ~] = findAverageAndStdevFiles(folderPath);
        field = readDavisField(avgFile, true);
        thisField = parameterDavisFieldToTable(field, pivCase, cfg, ...
            parameterGroup, parameterLabel, nSnapshots, deltaTLabel);
        parameterFields = [parameterFields; thisField]; %#ok<AGROW>

        valid = thisField.valid;
        summaryRow = table( ...
            string(parameterGroup), string(parameterLabel), string(pivCase.case_id), ...
            pivCase.aoa_deg, nSnapshots, string(deltaTLabel), nnz(valid), ...
            height(thisField), nnz(valid) / height(thisField), ...
            mean(thisField.u_streamwise_ms(valid), 'omitnan'), ...
            std(thisField.u_streamwise_ms(valid), 'omitnan'), ...
            mean(thisField.speed_ms(valid), 'omitnan'), ...
            'VariableNames', {'parameter_group', 'parameter_label', 'case_id', ...
            'aoa_deg', 'n_snapshots', 'delta_t_label', 'n_valid_vectors', ...
            'n_total_vectors', 'valid_fraction', 'mean_u_streamwise_ms', ...
            'spatial_std_u_streamwise_ms', 'mean_speed_ms'});
        parameterSummary = [parameterSummary; summaryRow]; %#ok<AGROW>
    end
end

function result = parameterDavisFieldToTable(field, pivCase, cfg, ...
        parameterGroup, parameterLabel, nSnapshots, deltaTLabel)
    u = cfg.pivStreamwiseSign * field.Vx;
    v = field.Vy;
    speed = hypot(u, v);
    result = table( ...
        repmat(string(parameterGroup), numel(field.X), 1), ...
        repmat(string(parameterLabel), numel(field.X), 1), ...
        repmat(string(pivCase.case_id), numel(field.X), 1), ...
        repmat(pivCase.aoa_deg, numel(field.X), 1), ...
        repmat(nSnapshots, numel(field.X), 1), ...
        repmat(string(deltaTLabel), numel(field.X), 1), ...
        field.X(:), field.Y(:), field.Y(:) - cfg.trailingEdgePosition_mm, ...
        u(:), v(:), speed(:), field.valid(:), ...
        'VariableNames', {'parameter_group', 'parameter_label', 'case_id', ...
        'aoa_deg', 'n_snapshots', 'delta_t_label', 'x_mm', 'y_mm', ...
        'y_relative_mm', 'u_streamwise_ms', 'v_normal_ms', 'speed_ms', 'valid'});
end

function [avgFile, stdevFile] = findAverageAndStdevFiles(folderPath)
    avgCandidates = {'avg.dat', 'Avg.dat', 'B00001.dat'};
    stdevCandidates = {'stdev.dat', 'RMS.dat', 'B00002.dat'};
    avgFile = findFirstExisting(folderPath, avgCandidates);
    stdevFile = findFirstExisting(folderPath, stdevCandidates);
    if avgFile == "" || stdevFile == ""
        error('Could not identify average and stdev files in %s.', folderPath);
    end
end

function path = findFirstExisting(folderPath, names)
    path = "";
    for iName = 1:numel(names)
        candidate = fullfile(folderPath, names{iName});
        if isfile(candidate)
            path = string(candidate);
            return;
        end
    end
end

function field = readDavisField(filePath, trimEdge)
    filePath = char(filePath);
    [nCols, nRows] = extractDavisGridSize(filePath);
    data = importdata(filePath, ' ', 3);
    matrix = data.data;

    field.X = reshape(matrix(:, 1), nCols, nRows).';
    field.Y = reshape(matrix(:, 2), nCols, nRows).';
    field.Vx = reshape(matrix(:, 3), nCols, nRows).';
    field.Vy = reshape(matrix(:, 4), nCols, nRows).';
    field.valid = reshape(matrix(:, 5), nCols, nRows).' > 0;

    if trimEdge
        field.X = field.X(1:end - 1, 1:end - 1);
        field.Y = field.Y(1:end - 1, 1:end - 1);
        field.Vx = field.Vx(1:end - 1, 1:end - 1);
        field.Vy = field.Vy(1:end - 1, 1:end - 1);
        field.valid = field.valid(1:end - 1, 1:end - 1);
    end
end

function [nCols, nRows] = extractDavisGridSize(filePath)
    fid = fopen(filePath, 'r');
    if fid < 0
        error('Cannot open %s.', filePath);
    end
    cleanup = onCleanup(@() fclose(fid));
    fgetl(fid);
    fgetl(fid);
    thirdLine = fgetl(fid);
    tokens = regexp(thirdLine, 'I=(\d+),\s*J=(\d+)', 'tokens', 'once');
    if isempty(tokens)
        error('Could not parse DaVis grid size from %s.', filePath);
    end
    nCols = str2double(tokens{1});
    nRows = str2double(tokens{2});
end

function result = davisFieldToTable(field, pivCase, processing, fieldType, cfg)
    u = cfg.pivStreamwiseSign * field.Vx;
    v = field.Vy;
    speed = hypot(u, v);
    result = table( ...
        repmat(string(pivCase.case_id), numel(field.X), 1), ...
        repmat(pivCase.aoa_deg, numel(field.X), 1), ...
        repmat(string(processing), numel(field.X), 1), ...
        repmat(string(fieldType), numel(field.X), 1), ...
        field.X(:), field.Y(:), field.Y(:) - cfg.trailingEdgePosition_mm, ...
        u(:), v(:), speed(:), field.valid(:), ...
        'VariableNames', {'case_id', 'aoa_deg', 'processing', 'field_type', ...
        'x_mm', 'y_mm', 'y_relative_mm', 'u_streamwise_ms', 'v_normal_ms', ...
        'speed_ms', 'valid'});
end

function result = davisRmsFieldToTable(field, pivCase, processing, cfg)
    uRms = abs(field.Vx);
    vRms = abs(field.Vy);
    result = table( ...
        repmat(string(pivCase.case_id), numel(field.X), 1), ...
        repmat(pivCase.aoa_deg, numel(field.X), 1), ...
        repmat(string(processing), numel(field.X), 1), ...
        field.X(:), field.Y(:), field.Y(:) - cfg.trailingEdgePosition_mm, ...
        uRms(:), vRms(:), hypot(uRms(:), vRms(:)), field.valid(:), ...
        'VariableNames', {'case_id', 'aoa_deg', 'processing', ...
        'x_mm', 'y_mm', 'y_relative_mm', 'u_rms_ms', 'v_rms_ms', ...
        'speed_rms_ms', 'valid'});
end

function profile = extractPivProfile(cfg, pivCase, avgField, stdevField, processing)
    xColumns = mean(avgField.X, 1, 'omitnan');
    [~, selectedColumn] = min(abs(xColumns - cfg.pivProfileX_mm));
    valid = avgField.valid(:, selectedColumn) & stdevField.valid(:, selectedColumn);

    y = avgField.Y(:, selectedColumn);
    yRelative = y - cfg.trailingEdgePosition_mm;
    u = cfg.pivStreamwiseSign * avgField.Vx(:, selectedColumn);
    uRms = abs(stdevField.Vx(:, selectedColumn));

    profile = table( ...
        repmat(string(pivCase.case_id), numel(y), 1), ...
        repmat(pivCase.aoa_deg, numel(y), 1), ...
        repmat(string(processing), numel(y), 1), ...
        repmat(xColumns(selectedColumn), numel(y), 1), ...
        y, yRelative, u, uRms, valid, ...
        'VariableNames', {'case_id', 'aoa_deg', 'processing', ...
        'x_profile_mm', 'y_mm', 'y_relative_mm', ...
        'mean_velocity_ms', 'rms_velocity_ms', 'valid'});
end

function summary = summarizePivProcessingFolder(cfg, pivCase)
    caseDir = fullfile(cfg.pivProcessedDir, pivCase.processedFolder);
    folders = dir(fullfile(caseDir, '*AvgStDev'));
    summary = table();
    for iFolder = 1:numel(folders)
        folderPath = fullfile(folders(iFolder).folder, folders(iFolder).name);
        try
            [avgFile, stdevFile] = findAverageAndStdevFiles(folderPath);
            avgField = readDavisField(avgFile, true);
            stdevField = readDavisField(stdevFile, true);
        catch
            continue;
        end

        valid = avgField.valid;
        u = cfg.pivStreamwiseSign * avgField.Vx;
        v = avgField.Vy;
        uRms = abs(stdevField.Vx);
        vRms = abs(stdevField.Vy);

        row = table( ...
            string(pivCase.case_id), pivCase.aoa_deg, ...
            erase(string(folders(iFolder).name), 'AvgStDev'), ...
            nnz(valid), numel(valid), nnz(valid) / numel(valid), ...
            mean(u(valid), 'omitnan'), mean(v(valid), 'omitnan'), ...
            mean(hypot(u(valid), v(valid)), 'omitnan'), ...
            mean(uRms(valid), 'omitnan'), mean(vRms(valid), 'omitnan'), ...
            'VariableNames', {'case_id', 'aoa_deg', 'processing', ...
            'n_valid_vectors', 'n_total_vectors', 'valid_fraction', ...
            'mean_u_streamwise_ms', 'mean_v_normal_ms', 'mean_speed_ms', ...
            'spatial_mean_u_rms_ms', 'spatial_mean_v_rms_ms'});
        summary = [summary; row]; %#ok<AGROW>
    end
end

function [selfFields, windowSummary, comparePoints, compareSummary] = processSelfPiv(cfg)
    pixelToMm = calibratePivImage(cfg.selfPiv.calibrationImagePath);
    selfFields = table();
    windowSummary = table();

    for windowSize = cfg.selfPiv.windowSizes_px
        field = computeSelfPivField(cfg.selfPiv.imagePath, windowSize, ...
            cfg.selfPiv.deltaT_s, pixelToMm, cfg.pivStreamwiseSign);
        selfFields = [selfFields; field]; %#ok<AGROW>

        valid = true(height(field), 1);
        row = table( ...
            windowSize, height(field), ...
            mean(field.u_streamwise_ms(valid), 'omitnan'), ...
            mean(field.v_normal_ms(valid), 'omitnan'), ...
            mean(field.speed_ms(valid), 'omitnan'), ...
            std(field.u_streamwise_ms(valid), 'omitnan'), ...
            'VariableNames', {'window_size_px', 'n_vectors', ...
            'mean_u_streamwise_ms', 'mean_v_normal_ms', 'mean_speed_ms', ...
            'spatial_std_u_streamwise_ms'});
        windowSummary = [windowSummary; row]; %#ok<AGROW>
    end

    reference = readDavisField(cfg.selfPiv.referenceDavisPath, true);
    self32 = selfFields(selfFields.window_size_px == 32, :);
    selfU = reshape(self32.u_streamwise_ms, size(reference.X));
    selfV = reshape(self32.v_normal_ms, size(reference.X));
    davisU = cfg.pivStreamwiseSign * reference.Vx;
    davisV = reference.Vy;
    valid = reference.valid;

    du = selfU - davisU;
    dv = selfV - davisV;
    comparePoints = table( ...
        reference.X(:), reference.Y(:), reference.Y(:) - cfg.trailingEdgePosition_mm, ...
        davisU(:), davisV(:), selfU(:), selfV(:), du(:), dv(:), valid(:), ...
        'VariableNames', {'x_mm', 'y_mm', 'y_relative_mm', ...
        'davis_u_streamwise_ms', 'davis_v_normal_ms', ...
        'self_u_streamwise_ms', 'self_v_normal_ms', ...
        'difference_u_ms', 'difference_v_ms', 'valid'});

    compareSummary = table( ...
        sqrt(mean(du(valid).^2, 'omitnan')), ...
        sqrt(mean(dv(valid).^2, 'omitnan')), ...
        mean(abs(du(valid)), 'omitnan'), ...
        mean(abs(dv(valid)), 'omitnan'), ...
        mean(davisU(valid), 'omitnan'), ...
        mean(selfU(valid), 'omitnan'), ...
        nnz(valid), ...
        'VariableNames', {'rmse_u_ms', 'rmse_v_ms', 'mean_abs_u_ms', ...
        'mean_abs_v_ms', 'mean_davis_u_ms', 'mean_self_u_ms', 'n_valid'});
end

function pixelToMm = calibratePivImage(calibrationImagePath)
    imread(calibrationImagePath);
    x1 = 330;
    x2 = 1050;
    y1 = 1160;
    y2 = 1150;
    knownDistance_mm = 80;
    pixelToMm = knownDistance_mm / hypot(x2 - x1, y2 - y1);
end

function field = computeSelfPivField(imagePath, windowSize, deltaT_s, pixelToMm, streamwiseSign)
    [image1, image2] = splitPivImage(imagePath);
    [height, width] = size(image1);
    nRows = floor(height / windowSize);
    nCols = floor(width / windowSize);
    nVectors = nRows * nCols;

    rowIndex = zeros(nVectors, 1);
    colIndex = zeros(nVectors, 1);
    x_mm = zeros(nVectors, 1);
    y_mm = zeros(nVectors, 1);
    vx = zeros(nVectors, 1);
    vy = zeros(nVectors, 1);

    k = 0;
    for iCol = 1:nCols
        for iRow = 1:nRows
            k = k + 1;
            xPixels = (iCol - 1) * windowSize + 1:iCol * windowSize;
            yPixels = (iRow - 1) * windowSize + 1:iRow * windowSize;
            window1 = double(image1(yPixels, xPixels));
            window2 = double(image2(yPixels, xPixels));
            window1 = window1 - mean(window1(:));
            window2 = window2 - mean(window2(:));

            correlationMap = xcorr2(window1, window2);
            [~, vectorizedIndex] = max(correlationMap(:));
            [peakY, peakX] = ind2sub(size(correlationMap), vectorizedIndex);
            dx_px = peakX - windowSize;
            dy_px = peakY - windowSize;

            rowIndex(k) = iRow;
            colIndex(k) = iCol;
            x_mm(k) = ((iCol - 0.5) * windowSize) * pixelToMm;
            y_mm(k) = ((nRows - iRow + 0.5) * windowSize) * pixelToMm;
            vx(k) = dx_px * pixelToMm / deltaT_s / 1000;
            vy(k) = dy_px * pixelToMm / deltaT_s / 1000;
        end
    end

    u = streamwiseSign * vx;
    field = table( ...
        repmat(windowSize, nVectors, 1), rowIndex, colIndex, x_mm, y_mm, ...
        vx, vy, u, vy, hypot(u, vy), ...
        'VariableNames', {'window_size_px', 'row_index', 'col_index', ...
        'x_mm', 'y_mm', 'vx_image_ms', 'vy_image_ms', ...
        'u_streamwise_ms', 'v_normal_ms', 'speed_ms'});
end

function [image1, image2] = splitPivImage(imagePath)
    image = imread(imagePath);
    halfHeight = floor(size(image, 1) / 2);
    image2 = image(1:halfHeight, :);
    image1 = image(halfHeight + 1:2 * halfHeight, :);
end

function comparison = buildTechniqueComparisonProfile(cfg, hwaProfiles, pivProfiles)
    hwa = table( ...
        repmat("HWA", height(hwaProfiles), 1), ...
        string(hwaProfiles.aoa_deg), hwaProfiles.aoa_deg, ...
        hwaProfiles.y_relative_mm, hwaProfiles.mean_velocity_ms, ...
        hwaProfiles.rms_velocity_ms, true(height(hwaProfiles), 1), ...
        'VariableNames', {'technique', 'case_id', 'aoa_deg', ...
        'y_relative_mm', 'mean_velocity_ms', 'rms_velocity_ms', 'valid'});

    pivValid = pivProfiles.valid;
    piv = table( ...
        repmat("PIV", height(pivProfiles), 1), ...
        pivProfiles.case_id, pivProfiles.aoa_deg, ...
        pivProfiles.y_relative_mm, pivProfiles.mean_velocity_ms, ...
        pivProfiles.rms_velocity_ms, pivValid, ...
        'VariableNames', {'technique', 'case_id', 'aoa_deg', ...
        'y_relative_mm', 'mean_velocity_ms', 'rms_velocity_ms', 'valid'});

    comparison = [hwa; piv];

    pressurePath = fullfile(cfg.dataDir, 'pressure_probe_profile.csv');
    if isfile(pressurePath)
        pressure = readtable(pressurePath);
        if all(ismember({'aoa_deg', 'y_relative_mm', 'mean_velocity_ms'}, pressure.Properties.VariableNames))
            pressureRows = table( ...
                repmat("Pressure probe", height(pressure), 1), ...
                string(pressure.aoa_deg), pressure.aoa_deg, pressure.y_relative_mm, ...
                pressure.mean_velocity_ms, nan(height(pressure), 1), ...
                true(height(pressure), 1), ...
                'VariableNames', comparison.Properties.VariableNames);
            comparison = [comparison; pressureRows]; %#ok<AGROW>
        end
    end
end

function missingInputs = detectMissingInputs(cfg)
    pressurePath = fullfile(cfg.dataDir, 'pressure_probe_profile.csv');
    missingInputs = table( ...
        ["pressure_probe_profile"], ...
        ["Required for the report item comparing PIV, HWA, and pressure probe mean profiles."], ...
        [~isfile(pressurePath)], ...
        [string(pressurePath)], ...
        'VariableNames', {'input_name', 'why_needed', 'is_missing', 'expected_path'});
end
