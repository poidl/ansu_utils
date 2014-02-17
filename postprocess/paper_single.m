% plotting script
clear all;
close all;

run1=570;
fname=['../exp',num2str(run1),'/data/iteration_history.mat'];
varname= 'pns_hist';
load(fname, varname);
load(['../exp',num2str(run1),'/data/input_data.mat'], 'lats','longs');

load(['../exp',num2str(run1),'/data/stationindex.mat']);

lat=squeeze(lats(1,:,1));
lon=squeeze(longs(1,1,:));

vv1=pns_hist; % variable to plot
it=size(vv1,1);
%it=5;

va=squeeze(vv1(it,:,:));

% if 1
%     load(['../exp',num2str(run1),'/data/stationindex.mat']);
%     addpath(genpath(['../exp',num2str(run1),'/ansu/']));
%     regions=find_regions(va);
%     for iregion=1:length(regions)
%         region=regions{iregion};     
%         if ~ismember(istation,region)
%             va(region)=nan;
%         end
%     end
% end



sz=1.0*[14 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 
h=imagesc(lon,lat,va);
set(gca,'YDir','normal')
set(h,'alphadata',~isnan(va))
%     cmp=colormap(hot(128));
%     cmp2=cmp(find(cmp(:,1)==1,1,'first')-15:end,:);
%     colormap([fliplr(cmp2);flipud(cmp2)]) ;
%     maxi=max(abs(va(:)));
%     caxis([-maxi maxi])
col=colorbar()
yl=get(col,'ylabel');
ylpos=get(yl,'position');
ylpos(1)=ylpos(1)+2;
%set(yl,'position',ylpos)
set(yl,'string','Pressure [db]','rotation',270,'position',ylpos)

hold on
[jis,iis]=ind2sub(size(va),istation);
disp(['istation: ',num2str(lon(iis)),' Lon ',num2str(lat(jis)),' Lat'])
disp(['value at istation: ',num2str(va(istation),20)])
plot(lon(iis),lat(jis),'kx','markersize',20,'linewidth',3)

xlabel('Longitude')
ylabel('Latitude')

load ('../external_scripts/coast/coast_data.mat');
plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);

print('-dpng','-r400',['figures/paper_single'])
