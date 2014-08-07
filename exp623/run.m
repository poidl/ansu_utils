% add library paths
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_04'))
addpath(genpath('../external_scripts'))
addpath(genpath('.'))

% subsample data
%subsample_gk_ak_gamma;
sinusoidal_single_surf

save_netcdf_input_single_surf(sns,ctns,pns);
disp('writing done')
[drhodx,drhody]=delta_tilde_rho(sns,ctns,pns);
regions={[1:length(sns(:))]'};
[drho,res]=solve_lsqr_test(regions, drhodx, drhody);

