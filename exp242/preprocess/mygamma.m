% plotting script
clear all;
close all;

fname='../../exp238/data/iteration_history.mat';
varname= 'pns_hist';
load(fname, varname);
load('../../exp238/data/input_data.mat', 'lats','longs');

lat=squeeze(lats(1,:,1));
lon=squeeze(longs(1,1,:));

rpsurf=squeeze(pns_hist(1,:,:)); % variable to plot

fname='../../exp239/data/iteration_history.mat';
varname= 'pns_hist';
load(fname, varname);

gam_top=squeeze(pns_hist(end,:,:)); 

fname='../../exp240/data/iteration_history.mat';
varname= 'pns_hist';
load(fname, varname);

gam_bot=squeeze(pns_hist(end,:,:));

gam_rp=(rpsurf-gam_top)./(gam_bot-gam_top);
gam=gam_rp;


save('../data/gamma.mat','gam')

        sz=1.0*[14 10];
        figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 
h=imagesc(lon,lat,gam_rp);
set(gca,'YDir','normal')
set(h,'alphadata',~isnan(gam_rp))
        cmp=colormap(hot(128));
        cmp2=cmp(find(cmp(:,1)==1,1,'first')-15:end,:);
    colormap([fliplr(cmp2);flipud(cmp2)]) ;
    maxi=max(abs(gam_rp(:)));
    caxis([-maxi maxi])
colorbar
hold on
load ('../../external_scripts/coast/coast_data.mat');
plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);
title('gamma on surface')
%print('-dpdf','-r400','figures/compare_pressure_circshift.pdf')
