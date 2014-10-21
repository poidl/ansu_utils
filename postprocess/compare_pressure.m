% plotting script
clear all;
close all;

run1=642;
run2=643;
fname=['../exp',num2str(run1),'/data/iteration_history.mat'];
varname= 'pns_hist';
load(fname, varname);
%load(['../exp',num2str(run1),'/data/input_data.mat'], 'lats','longs');
gamma_i_exp=24;
load(['/home/nfs/z3439823/mymatlab/gamma_i/gamma_i_utils/exp0',num2str(gamma_i_exp),'/data/input_data.mat'])
lat=squeeze(lats(1,:,1));
lon=squeeze(longs(1,1,:));

vv1=pns_hist; % variable to plot

fname=['../exp',num2str(run2),'/data/iteration_history.mat'];
varname= 'pns_hist';
load(fname, varname);

vv2=pns_hist; % variable to plot

vv1it=squeeze(vv1(end,:,:));
if 0
    va=vv1it;
    h=imagesc(va)
    set(h,'alphadata',~isnan(va))
    set(gca,'YDir','normal')
    colorbar()
    keyboard
end
vv2it=squeeze(vv2(end,:,:));

vv2it(isnan(vv1it))=nan;
vv1it(isnan(vv2it))=nan;

diff=vv1it-vv2it;

%diff(abs(diff)>1e-8)=0;

if 1
    load(['../exp',num2str(run2),'/data/stationindex.mat']);
    addpath(genpath(['../exp',num2str(run1),'/ansu/']));
    regions=find_regions(diff);
    for iregion=1:length(regions)
        region=regions{iregion};     
        if ~ismember(istation,region)
            diff(region)=nan;
        end
    end
end

logp=true;

sz=1.0*[14 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 
if ~logp
    h=imagesc(lon,lat,diff);
else
    h=imagesc(lon,lat,log10(abs(diff)));
end
set(gca,'YDir','normal')
set(h,'alphadata',~isnan(diff))
if ~logp
    cmp=colormap(hot(128));
    cmp2=cmp(find(cmp(:,1)==1,1,'first')-15:end,:);
    colormap([fliplr(cmp2);flipud(cmp2)]) ;
    maxi=max(abs(diff(:)));
    caxis([-maxi maxi])
    colorbar
else
    colorbar
end

hold on
[jis,iis]=ind2sub(size(diff),istation)
% disp(['value at istation: ',num2str(diff(istation))])
% disp(['value at istation: ',num2str(vv1it(istation),20)])
plot(lon(iis),lat(jis),'b*','markersize',30,'linewidth',3)
plot(lon(iis),lat(jis),'ro','markersize',30,'linewidth',3)
plot(lon(iis),lat(jis),'go','markersize',40,'linewidth',3)

ilat1=18;
ilon1=48;
ilat2=5;
ilon2=85;

disp(['dp at bb 1: ',num2str(diff(ilat1,ilon1))])
disp(['dp at bb 2: ',num2str(diff(ilat2,ilon2))])
disp(['max(dp): ',num2str(max(abs(diff(:))))])
%keyboard
load ('../external_scripts/coast/coast_data.mat');
plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);
if ~logp
    title('Pressure difference [dbar]')
else
    title(['Pressure difference [log10(abs(dbar))]. Initial: ',num2str(vv2it(istation)),' dbar'])
end
print('-dpdf','-r400',['figures/compare_pressure_finals_',num2str(run1),'_',num2str(run2),'.pdf'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

