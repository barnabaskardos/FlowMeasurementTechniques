% PIV.m
% Compatibility entry point for the report-ready processing pipeline.
%
% The self-made PIV implementation, DaVis comparison, PIV statistics, and
% HWA processing are implemented in analysis/process_data.m. Run this file
% from MATLAB to regenerate the CSV tables under outputs/tables.

rootDir = fileparts(mfilename('fullpath'));
run(fullfile(rootDir, 'analysis', 'process_data.m'));
