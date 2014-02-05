% plotting script
clear all;
close all;

run1=568;
run2=568;
fname=['../exp',num2str(run1),'/data/iteration_history.mat'];
varname= 'pns_hist';
load(fname, varname);
load(['../exp',num2str(run1),'/data/input_data.mat'], 'lats','longs');

lat=squeeze(lats(1,:,1));
lon=squeeze(longs(1,1,:));

vv1=pns_hist; % variable to plot
it=size(vv1,1);
%it=5;
fname=['../exp',num2str(run2),'/data/iteration_history.mat'];
varname= 'pns_hist';
load(fname, varname);

vv2=pns_hist; % variable to plot

vv1it=squeeze(vv1(it,:,:));
vv2it=squeeze(vv2(it,:,:));

vv2it(isnan(vv1it))=nan;
vv1it(isnan(vv2it))=nan;

diff=vv1it-vv2it;

%diff(abs(diff)>1e-8)=0;

if 1
    load(['../exp',num2str(run1),'/data/stationindex.mat']);
    addpath(genpath(['../exp',num2str(run1),'/ansu/']));
    regions=find_regions(diff);
    for iregion=1:length(regions)
        region=regions{iregion};     
        if ~ismember(istation,region)
            diff(region)=nan;
        end
    end
end



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
[jis,iis]=ind2sub(size(diff),istation);
disp(['value at istation: ',num2str(diff(istation))])
disp(['value at istation: ',num2str(vv1it(istation),20)])
plot(lon(iis),lat(jis),'rx','markersize',20,'linewidth',2)

load ('../external_scripts/coast/coast_data.mat');
plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);
title('Pressure difference [db]')
print('-dpdf','-r400',['figures/compare_pressure_finals_',num2str(run1),'_',num2str(run2),'.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
final=repmat(vv1(end,:,:),[size(vv1,1),1,1]);
diff=abs(final-vv2);
diff=max(diff,[],3);
diff=max(diff,[],2);

xax=[0:size(vv1,1)-1];

figure()
semilogy(xax,diff)
hold on 
semilogy(xax,diff,'o')
print('-dpdf','-r400',['figures/compare_pressure_max_',num2str(run1),'_',num2str(run2),'.pdf'])
