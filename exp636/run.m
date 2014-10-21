% add library paths
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_04'))
addpath(genpath('../external_scripts'))
addpath(genpath('.'))

% subsample data
%subsample_gk_ak_gamma;
%sinusoidal

% construct zero-helicity ocean from the subsampled data
%zero_hel_no_continents

% run Paul's script
%omega_ak

postpr_gi

%plot
delete figures/*pdf
delete figures/*png
pressure_ns
%delta_pressure_ns
%drho
% if 0;
%     error_measure
%     error_measure_phi
% else
%     drho
%     error_measure_drho
     error_measure_derr_derrxy
% end
% 
% get_mean_pressure

