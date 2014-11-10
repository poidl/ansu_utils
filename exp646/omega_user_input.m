%datapath='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp644/data/';
datapath='./data/';
zonally_periodic=true;
clamp_on_point=true; 
error_measure='drho_local';
%error_measure='slope_difference'; % less stable than 'drho_local', needs damping. Need to monitor convergence, to be sure increase nit_after_wetting to 60.
no_land_mask=true; % 'true' may only work with the original Jackatt/McDougall 96 data set
use_b=false; % 'true' may not work.
solver = 'iterative';
delta=1e-11;
%weights='drho_local_n2'; % only works if error_measure='drho_local';
weights='';
    

if strcmp(error_measure,'drho_local')
    nit_after_wetting=15;
elseif strcmp(error_measure,'slope_difference')
    nit_after_wetting=60;
end
% number of iterations after the first iteration at which number of added points is zero.
%nit_max=30; % in case wetting keeps adding points for some reason, set a maximum number of iterations (will raise error)

save_iterations=true; % save variables for postprocessing 
history_file=[datapath,'iteration_history.mat']; 
