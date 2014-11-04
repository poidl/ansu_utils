%datapath='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp644/data/';
datapath='./data/';
use_b=false;
no_land_mask=false;
clamp_on_point=true;
zonally_periodic=true;
solver = 'iterative';
delta=1e-11;
%error_measure='slope_difference';
error_measure='drho_local';
    %weights='drho_local_n2'; % only works if error_measure='drho_local';
    weights='';
    
% nit = 6*1-1;
% nit = 6*1-1;
nit_after_wetting=15; % number of iterations after the first iteration at which number of added points is zero.
%nit_max=30; % in case wetting keeps adding points for some reason, set a maximum number of iterations (will raise error)

save_iterations=true; % save variables for postprocessing 
history_file=[datapath,'iteration_history.mat']; 
