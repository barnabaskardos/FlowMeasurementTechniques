
clear;
clc;
close all;

baseDir = fileparts(mfilename('fullpath'));
hwaDir = fullfile(baseDir, 'data', 'HWA', 'Group15');
calibrationPath = fullfile(hwaDir, 'calibration.txt');

calibrationData = readmatrix(calibrationPath);


calibrationData = calibrationData(:, 1:2);
calibrationData = calibrationData(all(isfinite(calibrationData), 2), :);

uCal = calibrationData(:, 1);
signalCal = calibrationData(:, 2);

% calibration.txt contains repeated signal levels; average them so interp1 is valid.
[signalUnique, ~, groupIdx] = unique(signalCal, 'sorted');
uUnique = accumarray(groupIdx, uCal, [], @mean);


dataFiles = dir(fullfile(hwaDir, 'aoa*_10.txt'));

records = struct('aoa', {}, 'position_mm', {}, 'mean_velocity', {}, ...
                 'rms_velocity', {}, 'sample_rate_hz', {}, ...
                 'nSamples', {}, 'fileName', {});

for iFile = 1:numel(dataFiles)
    fileName = dataFiles(iFile).name;
    tokens = regexp(fileName, '^aoa(-?\d+)_([0-9]+)_10\.txt$', 'tokens', 'once');
    if isempty(tokens)
        continue;
    end

    aoa = str2double(tokens{1});
    position_mm = str2double(tokens{2});
    filePath = fullfile(dataFiles(iFile).folder, fileName);

    rawData = readmatrix(filePath);
   

    probeSignal = rawData(:, 2);
    probeSignal = probeSignal(isfinite(probeSignal));


    velocity = interp1(signalUnique, uUnique, probeSignal, 'pchip', 'extrap');
    meanVelocity = mean(velocity, 'omitnan');
    velocityFluctuation = velocity - meanVelocity;
    rmsVelocity = sqrt(mean(velocityFluctuation.^2, 'omitnan'));

    time = rawData(:, 1);
    sampleRate_hz = 1 / mean(diff(time(isfinite(time))));

    records(end + 1).aoa = aoa; %#ok<SAGROW>
    records(end).position_mm = position_mm;
    records(end).mean_velocity = meanVelocity;
    records(end).rms_velocity = rmsVelocity;
    records(end).sample_rate_hz = sampleRate_hz;
    records(end).nSamples = numel(velocity);
    records(end).fileName = fileName;
end

if isempty(records)
    error('No valid HWA records could be processed.');
end

T = struct2table(records);
T = sortrows(T, {'aoa', 'position_mm'});

aoaValues = unique(T.aoa);
commonPositions = unique(T.position_mm);
for iAoa = 1:numel(aoaValues)
    posThisAoa = unique(T.position_mm(T.aoa == aoaValues(iAoa)));
    commonPositions = intersect(commonPositions, posThisAoa);
end

if isempty(commonPositions)
    error('No common measurement positions found across AoA datasets.');
end

Tplot = T(ismember(T.position_mm, commonPositions), :);
Tplot = sortrows(Tplot, {'position_mm', 'aoa'});

figure('Color', 'w');
hold on;
grid on;
box on;

lineColors = lines(numel(aoaValues));

for iAoa = 1:numel(aoaValues)
    thisAoa = aoaValues(iAoa);
    idx = Tplot.aoa == thisAoa;
    x = Tplot.position_mm(idx);
    y = Tplot.mean_velocity(idx);

    [xSorted, order] = sort(x);
    ySorted = y(order);

    plot(xSorted, ySorted, '-o', ...
        'Color', lineColors(iAoa, :), ...
        'LineWidth', 1.6, ...
        'MarkerSize', 5, ...
        'DisplayName', sprintf('AoA %d deg', thisAoa));
end

xlabel('Position [mm]');
ylabel('Mean Velocity [m/s]');
title('HWA Mean Velocity vs Position (Interpolated from 22-point calibration)');
legend('Location', 'best');
xline(55, '--k', '55 mm', 'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'center');

hold off;

spectraWindowSamples = 16384;
spectraOverlapFraction = 0.5;
spectraPlotMinFrequency_hz = 1;
spectraPlotMaxFrequency_hz = 5000;
minimumPeakFrequency_hz = 10;
minimumPeakSeparation_hz = 20;
nDominantPeaks = 3;

spectraAoAValues = [0, 5, 15];
trailingEdgePosition_mm = 55;
wakeSpectraPositions_mm = [trailingEdgePosition_mm, ...
    trailingEdgePosition_mm + 4, trailingEdgePosition_mm - 4];
wakeSpectraLabels = {'Trailing edge', 'Above TE', 'Below TE'};

spectraSummary = table('Size', [0 8], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double', ...
                      'double', 'double', 'double'}, ...
    'VariableNames', {'Region', 'Position_mm', 'AoA_deg', 'RMSVelocity', ...
                      'Peak1_Hz', 'Peak2_Hz', 'Peak3_Hz', ...
                      'FrequencyResolution_Hz'});

spectraColors = lines(numel(spectraAoAValues));

for iPosition = 1:numel(wakeSpectraPositions_mm)
    thisPosition_mm = wakeSpectraPositions_mm(iPosition);
    thisRegion = wakeSpectraLabels{iPosition};

    figure('Color', 'w');
    hold on;
    grid on;
    box on;

    for iAoa = 1:numel(spectraAoAValues)
        thisAoa = spectraAoAValues(iAoa);
        filePath = fullfile(hwaDir, sprintf('aoa%d_%d_10.txt', ...
            thisAoa, thisPosition_mm));


        rawData = readmatrix(filePath);
        time = rawData(:, 1);
        probeSignal = rawData(:, 2);
        velocity = interp1(signalUnique, uUnique, probeSignal, 'pchip', 'extrap');
        validSamples = isfinite(time) & isfinite(velocity);
        time = time(validSamples);
        velocity = velocity(validSamples);

        [frequency_hz, phi_uu, frequencyResolution_hz] = computeEnergySpectrum( ...
            time, velocity, spectraWindowSamples, spectraOverlapFraction);

        plotMask = frequency_hz >= spectraPlotMinFrequency_hz & ...
            frequency_hz <= spectraPlotMaxFrequency_hz;
        loglog(frequency_hz(plotMask), phi_uu(plotMask), ...
            'Color', spectraColors(iAoa, :), ...
            'LineWidth', 1.4, ...
            'DisplayName', sprintf('AoA %d deg', thisAoa));

        peakFrequencies_hz = dominantFrequencies(frequency_hz, phi_uu, ...
            minimumPeakFrequency_hz, spectraPlotMaxFrequency_hz, ...
            minimumPeakSeparation_hz, nDominantPeaks);
        peakFrequencies_hz(end + 1:nDominantPeaks) = NaN;

        meanVelocity = mean(velocity, 'omitnan');
        rmsVelocity = sqrt(mean((velocity - meanVelocity).^2, 'omitnan'));

        spectraSummary(end + 1, :) = {string(thisRegion), thisPosition_mm, ...
            thisAoa, rmsVelocity, peakFrequencies_hz(1), ...
            peakFrequencies_hz(2), peakFrequencies_hz(3), ...
            frequencyResolution_hz}; %#ok<SAGROW>
    end

    xlabel('Frequency [Hz]');
    ylabel('\phi_{uu} [arb. units]');
    title(sprintf('HWA Energy Spectra at %s (y = %d mm)', ...
        thisRegion, thisPosition_mm));
    legend('Location', 'best');
    xlim([spectraPlotMinFrequency_hz, spectraPlotMaxFrequency_hz]);
    set(gca, 'YScale', 'log');
    hold off;
end

disp(spectraSummary);

fprintf('\nFFT/Welch window size effect for fs = %.1f Hz:\n', median(T.sample_rate_hz));
windowSamplesToCompare = [4096, 8192, 16384, 32768];
for iWin = 1:numel(windowSamplesToCompare)
    nWin = windowSamplesToCompare(iWin);
    fprintf('  Nwin = %5d samples: T = %.3f s, df = %.3f Hz, Nyquist = %.1f Hz\n', ...
        nWin, nWin / median(T.sample_rate_hz), ...
        median(T.sample_rate_hz) / nWin, median(T.sample_rate_hz) / 2);
end

function [frequency_hz, phi_uu, frequencyResolution_hz] = computeEnergySpectrum( ...
        time, velocity, windowSamples, overlapFraction)
    time = time(:);
    velocity = velocity(:);

    valid = isfinite(time) & isfinite(velocity);
    time = time(valid);
    velocity = velocity(valid);
    velocity = velocity - mean(velocity, 'omitnan');

    sampleRate_hz = 1 / mean(diff(time));
    windowSamples = min(windowSamples, numel(velocity));
    overlapSamples = floor(overlapFraction * windowSamples);
    nfft = windowSamples;

    window = hann(windowSamples, 'periodic');
    [phi_uu, frequency_hz] = pwelch(velocity, window, overlapSamples, ...
        nfft, sampleRate_hz);

    frequencyResolution_hz = sampleRate_hz / nfft;
end

function peakFrequencies_hz = dominantFrequencies(frequency_hz, phi_uu, ...
        minFrequency_hz, maxFrequency_hz, minSeparation_hz, nPeaks)
    valid = frequency_hz >= minFrequency_hz & frequency_hz <= maxFrequency_hz;
    candidate = valid;
    candidate(2:end - 1) = candidate(2:end - 1) & ...
        phi_uu(2:end - 1) > phi_uu(1:end - 2) & ...
        phi_uu(2:end - 1) >= phi_uu(3:end);
    candidate(1) = false;
    candidate(end) = false;

    peakFrequencies = frequency_hz(candidate);
    peakEnergies = phi_uu(candidate);
    [~, order] = sort(peakEnergies, 'descend');
    peakFrequencies = peakFrequencies(order);

    peakFrequencies_hz = [];
    for iPeak = 1:numel(peakFrequencies)
        farEnough = isempty(peakFrequencies_hz) || ...
            all(abs(peakFrequencies(iPeak) - peakFrequencies_hz) >= minSeparation_hz);
        if farEnough
            peakFrequencies_hz(end + 1) = peakFrequencies(iPeak); %#ok<AGROW>
        end

        if numel(peakFrequencies_hz) == nPeaks
            break;
        end
    end
end
