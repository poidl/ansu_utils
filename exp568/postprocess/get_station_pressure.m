clear all;

fname='data/iteration_history.mat';
varname= 'pns_hist';
load(fname, varname);

vv=pns_hist; % variable to plot
nit=size(vv,1);

load('data/stationindex.mat')
vv=vv(:,:);

for ii=1:nit;
    v=vv(ii,istation);
    disp(['surf. pressure at station ', num2str(v,16),'  it: ',num2str(ii-1)])
end