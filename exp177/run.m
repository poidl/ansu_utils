 % add library paths
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../external_scripts'))
addpath(genpath('.'))

% construct idealized data set
idealized_01

% construct zero-helicity ocean from the subsampled data
%zero_hel

% run Paul's script
%dbstop in omega_ak at 53
omega_ak

% plot
pressure_ns
delta_pressure_ns
%phi_prime
error_measure
error_measure_phi
get_mean_pressure
%corr_rIdeal_multiple_pics
%corr_rIdeal_phiprime_forLargePhi
