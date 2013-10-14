% plot scalar as function of iterations, e.g. rms slope error
clear all;

fname='data/iteration_history.mat';
varname= 'pns_hist';
load(fname, varname);

vv=pns_hist; % variable to plot
v1=vv(1,:,:);
v2=vv(end,:,:);
mp_initial=nanmean(v1(:));
mp_final=nanmean(v2(:));

disp(['mean pressure initial: ', num2str(mp_initial)])
disp(['mean pressure final: ', num2str(mp_final)])
