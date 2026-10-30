% Run this script from any MATLAB working directory. No extra toolbox needed.
folder=fileparts(mfilename('fullpath'));
addpath(folder);
verify_narrowband_validation(fullfile(folder,'narrowband_validation.mat'));
verify_narrowband_validation(fullfile(folder,'narrowband_boundary_valid.mat'));
verify_narrowband_validation(fullfile(folder,'narrowband_boundary_alias.mat'));
verify_continuous_ensemble(fullfile(folder,'continuous_ensemble_validation.mat'));
