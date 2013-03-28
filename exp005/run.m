% subsample 4x4
% 
clear all
close all
if 1
    % subsample data
    subsample_gk_ak_gamma;
    % construct zero-helicity ocean from the subsampled data
    zero_hel
    % run Paul's script
    gamma_ak
end
% plot
ansu_plot
error_measure



