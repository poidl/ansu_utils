use_b=true;
no_land_mask=true;
zonally_periodic=true;
solver = 'iterative';
delta=1e-11;
% nit = 6*1-1;
% nit = 6*1-1;
nit_after_wetting=3; % number of iterations after the first iteration at which number of added points is zero.
nit_max=6*3-1; % in case wetting keeps adding points for some reason, set a maximum number of iterations (will raise error)

save_iterations=true; % save variables for postprocessing 
history_file='data/iteration_history.mat'; % absolute path or path relative to the topmost calling function 
