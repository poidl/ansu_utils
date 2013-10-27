% plotting script
clear all;
close all;

fname='../exp238/data/iteration_history.mat';
varname= 'pns_hist';
load(fname, varname);
load('../exp238/data/input_data.mat', 'lats','longs');

lat=squeeze(lats(1,:,1));
lon=squeeze(longs(1,1,:));

vv1=pns_hist; % variable to plot

fname='../exp240/data/iteration_history.mat';
varname= 'pns_hist';
load(fname, varname);

vv2=pns_hist; % variable to plot

vv1=squeeze(vv1(1,:,:));
vv2=squeeze(vv2(end,:,:));

vv2(isnan(vv1))=nan;
vv1(isnan(vv2))=nan;

diff=vv1-vv2;

disp('any(diff(:)><0): ')
any(diff(:)<0)

        sz=1.0*[14 10];
        figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 
h=imagesc(lon,lat,diff);
set(gca,'YDir','normal')
set(h,'alphadata',~isnan(diff))
        cmp=colormap(hot(128));
        cmp2=cmp(find(cmp(:,1)==1,1,'first')-15:end,:);
    colormap([fliplr(cmp2);flipud(cmp2)]) ;
    maxi=max(abs(diff(:)));
    caxis([-maxi maxi])
colorbar
hold on
load ('../external_scripts/coast/coast_data.mat');
plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);
title('Pressure difference [db]')
print('-dpdf','-r400','figures/compare_pressure_circshift.pdf')

