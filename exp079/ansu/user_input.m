zonally_periodic=true;
choice= 'epsilon'; % density gradient error
solver = 'iterative';
nit = 21;

save_iterations=true; % save variables for postprocessing 
history_file='data/iteration_history.mat'; % absolute path or path relative to the topmost calling function 


glevels=[27.5]; % TODO: should be defined somewhere else