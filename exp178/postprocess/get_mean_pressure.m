% plot scalar as function of iterations, e.g. rms slope error
clear all;

fname='data/iteration_history.mat';
varname= 'pns_hist';
load(fname, varname);

vv=pns_hist(1,:,:); % variable to plot

mpressure=nanmean(vv(:));
disp(['mean pressure: ', num2str(mpressure)])


