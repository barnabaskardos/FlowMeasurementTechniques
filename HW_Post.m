% HW_Post.m
% Compatibility entry point for the report-ready processing pipeline.
%
% HWA calibration, sampling-time convergence, wake profiles, spectra, and
% the PIV comparison tables are implemented in analysis/process_data.m. Run
% this file from MATLAB to regenerate the CSV tables under outputs/tables.

rootDir = fileparts(mfilename('fullpath'));
run(fullfile(rootDir, 'analysis', 'process_data.m'));
