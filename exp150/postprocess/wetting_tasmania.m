% plotting script
% find the reason why wetting stopps around Tasmania in exp101
clear all;
close all;

restoredefaultpath
addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('../../external_scripts'))
addpath(genpath('../ansu'))
addpath(genpath('.'))


fname='../data/iteration_history.mat';
varname= 'sns_hist';
load(fname, varname);
varname= 'ctns_hist';
load(fname, varname);
varname= 'pns_hist';
load(fname, varname);

load('../data/input_data.mat', 'lats','longs');

lat=squeeze(lats(1,:,1));
lon=squeeze(longs(1,1,:));

vv1=pns_hist; % variable to plot
nit=size(vv1,1);

fname='../../exp099/data/iteration_history.mat';
varname= 'pns_hist';
load(fname, varname);

fname='../../exp099/data/input_data.mat';
varname= 's';
load(fname, varname);
varname= 'ct';
load(fname, varname);
varname= 'p';
load(fname, varname);


vv2=pns_hist; % variable to plot

itind=33;
pns=squeeze(vv1(itind,:,:));
sns=squeeze(sns_hist(itind,:,:));
ctns=squeeze(ctns_hist(itind,:,:));


jj=find(lat==-45);
ii=find(lon==144);

[glon,glat]=meshgrid(lon,lat);

en=jj+(ii-2)*size(s,2);
neighbour=jj+(ii-1)*size(s,2);

[a,b,c] = depth_ntp_iter(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,en),ct(:,en),p(:,en))


va=pns;

sz=1.0*[14 10];
        figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 
h=imagesc(lon,lat,va);
set(gca,'YDir','normal')
set(h,'alphadata',~isnan(va))
colorbar
hold on
load ('../../external_scripts/coast/coast_data.mat');
plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);
title('Pressure difference [db]')
print('-dpdf','-r400','../figures/compare_pressure.pdf')