% add library paths
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../external_scripts'))
addpath(genpath('.'))

% subsample data
subsample_gk_ak_gamma;

% construct zero-helicity ocean from the subsampled data
%zero_hel_no_continents

% run Paul's script
omega_ak

%plot
pressure_ns
%delta_pressure_final
% drho
% if 0;
%     error_measure
%     error_measure_phi
% else
     drho
%     error_measure_drho
     error_measure_drho_drhoxy
% end
% 
% get_mean_pressure

