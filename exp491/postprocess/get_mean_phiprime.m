% plot scalar as function of iterations, e.g. rms slope error
clear all;

fname='data/iteration_history.mat';
varname= 'phiprime_e_hist';
load(fname, varname);

vv=phiprime_e_hist; % variable to plot
nit=size(vv,1);


for ii=1:nit;
    v=vv(ii,:,:);
    mp=nanmean(v(:));
    disp(['mean phiprime ',num2str(ii-1), ' ', num2str(mp)])
end
