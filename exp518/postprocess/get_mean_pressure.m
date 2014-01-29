% plot scalar as function of iterations, e.g. rms slope error
clear all;

fname='data/iteration_history.mat';
varname= 'pns_hist';
load(fname, varname);

vv=pns_hist; % variable to plot
nit=size(vv,1);


for ii=1:nit;
    v=vv(ii,:,:);
    mp=nanmean(v(:));
    disp(['mean pressure ',num2str(ii-1), ' ', num2str(mp)])
end
