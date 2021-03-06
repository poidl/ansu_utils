 % add library paths
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../external_scripts'))
addpath(genpath('.'))

% subsample data
subsample_gk_ak_gamma;

% construct zero-helicity ocean from the subsampled data
zero_hel

% run Paul's script
gamma_ak

% plot
ansu_plot
error_measure



