no_land_mask=false;
zonally_periodic=true;
solver = 'iterative';
delta=1e-11;
nit = 6*1-1;
nit = 6*1-1;
nit=20; % should stop earlier
nit_after_wetting=3; % number of iterations after the first iteration at which number of added points is zero.

save_iterations=true; % save variables for postprocessing 
history_file='data/iteration_history.mat'; % absolute path or path relative to the topmost calling function 
