% plotting script
clear all;
close all;

restoredefaultpath
addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('../../external_scripts'))
addpath(genpath('.'))

fname='../../exp239/data/iteration_history.mat';
varname= {'sns_hist','ctns_hist','pns_hist'};
load(fname, varname{:});
load('../../exp239/data/input_data.mat', 'lats','longs');

lat=squeeze(lats(1,:,1));
lon=squeeze(longs(1,1,:));

s1=squeeze(sns_hist(end,:,:)); 
ct1=squeeze(ctns_hist(end,:,:));
p1=squeeze(pns_hist(end,:,:)); 

fname='../../exp240/data/iteration_history.mat';
varname= {'sns_hist','ctns_hist','pns_hist'};
load(fname, varname{:});

s2=squeeze(sns_hist(end,:,:)); 
ct2=squeeze(ctns_hist(end,:,:));
p2=squeeze(pns_hist(end,:,:)); 

pmid=0.5*(p1+p2);
b=1./(gsw_rho(s1,ct1,pmid)-gsw_rho(s2,ct2,pmid));

save('../data/b.mat','b')

v2p=b;
        sz=1.0*[14 10];
        figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 
h=imagesc(lon,lat,v2p);
set(gca,'YDir','normal')
set(h,'alphadata',~isnan(v2p))
      %  cmp=colormap(hot(128));
     %   cmp2=cmp(find(cmp(:,1)==1,1,'first')-15:end,:);
    %colormap([fliplr(cmp2);flipud(cmp2)]) ;
    %maxi=max(abs(v2p(:)));
    %caxis([-maxi maxi])
colorbar
hold on
load ('../../external_scripts/coast/coast_data.mat');
plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);
title('b')
print('-dpdf','-r400','../figures/myb.pdf')
