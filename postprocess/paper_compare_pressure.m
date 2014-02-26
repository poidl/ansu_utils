% plotting script
clear all;
close all;

run1=571;
run2=572;
fname=['../exp',num2str(run1),'/data/iteration_history.mat'];
varname= 'pns_hist';
load(fname, varname);
load(['../exp',num2str(run1),'/data/input_data.mat'], 'lats','longs');

lat=squeeze(lats(1,:,1));
lon=squeeze(longs(1,1,:));

vv1=pns_hist; % variable to plot
it=size(vv1,1);

fname=['../exp',num2str(run2),'/data/iteration_history.mat'];
varname= 'pns_hist';
load(fname, varname);

vv2=pns_hist; % variable to plot

vv1it=squeeze(vv1(it,:,:));
vv2it=squeeze(vv2(it,:,:));
diff=vv1it-vv2it;

load(['../exp',num2str(run1),'/data/stationindex.mat']);

sz=0.9*[15 18];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)])


subplot(2,1,1)
h1=imagesc(lon,lat,diff);
set(gca,'YDir','normal')
set(h1,'alphadata',~isnan(diff))
cb1=colorbar()
hold on
[jis,iis]=ind2sub(size(diff),istation);
disp(['value at istation: ',num2str(diff(istation))])
disp(['value at istation: ',num2str(vv1it(istation),20)])
plot(lon(iis),lat(jis),'kx','markersize',20,'linewidth',3)

xlabel('Longitude')
ylabel('Latitude')
load ('../external_scripts/coast/coast_data.mat');
plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(2,1,2)
h2=imagesc(lon,lat,vv1it);
set(gca,'YDir','normal')
set(h2,'alphadata',~isnan(vv1it))
cb2=colorbar()
hold on
[jis,iis]=ind2sub(size(vv1it),istation);
disp(['istation: ',num2str(lon(iis)),' Lon ',num2str(lat(jis)),' Lat'])
disp(['value at istation: ',num2str(vv1it(istation),20)])
plot(lon(iis),lat(jis),'kx','markersize',20,'linewidth',3)

xlabel('Longitude')
ylabel('Latitude')
load ('../external_scripts/coast/coast_data.mat');
plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmp=colormap(hot(128));
cmp2=cmp(find(cmp(:,1)==1,1,'first')-15:end,:);
cm1=colormap([fliplr(cmp2);flipud(cmp2)]) ;
cm2=jet;
cmap=[cm1;cm2];

colormap(cmap);

m=size(cm1,1);
va=diff;
%cmax1=max(abs(va(:)));
%cmin1=-cmax1;
cmax1=0.5;
cmin1=-cmax1;
g=min(m,round((m-1)*(0-cmin1)/(cmax1-cmin1))+1);
C1 = min(m, round( (m-1)*(va-cmin1)/(cmax1-cmin1)  )   +1);

va=vv1it;
cmin=min(va(:));
cmax=max(va(:));
cmin=0;
cmax=1800;
C2 = min(m,round((m-1)*(va-cmin)/(cmax-cmin))+1);
C2 = m+C2;

set(h1,'CData',C1);
set(h2,'CData',C2);

ax = findobj(gcf,'Type','axes');
set(ax,'CLim', [1 2*m])

set(cb1,'ylim',[1 m])
set(cb2,'ylim',[m+1 2*m])

nticks=5;
set(cb1,'ytick',linspace(1,m,nticks))
set(cb1,'yticklabel',linspace(cmin1,cmax1,nticks))

set(cb2,'ytick',linspace(m+1,2*m,nticks))
set(cb2,'yticklabel',linspace(cmin,cmax,nticks))


yl=get(cb1,'ylabel');
ylpos=get(yl,'position');
ylpos(1)=ylpos(1)+2;
%set(yl,'position',ylpos)
set(yl,'string','Pressure difference EXP1''-EXP5 [db]','rotation',270,'position',ylpos)


yl=get(cb2,'ylabel');
ylpos=get(yl,'position');
ylpos(1)=ylpos(1)+2;
%set(yl,'position',ylpos)
set(yl,'string','Pressure EXP1'' [db]','rotation',270,'position',ylpos)


text(-0.15,1.1,'b)','units','normalized','fontsize',15,'Parent',ax(2))
text(-0.15,1.1,'a)','units','normalized','fontsize',15,'Parent',ax(4))

%tick=get(cb1,'ytick')
%ltick=length(tick);
%set(cb1,'ytick',linspace(cmin1,cmax1,4))


%caxis([min(C1(:)) max(C2(:))])



% maxi=max(abs(diff(:)));
% caxis([-maxi maxi])
% col=colorbar()
% yl=get(col,'ylabel');
% ylpos=get(yl,'position');
% ylpos(1)=ylpos(1)+2;
% %set(yl,'position',ylpos)
% set(yl,'string','Pressure difference [db]','rotation',270,'position',ylpos)
% 
% 
% col=colorbar()
% yl=get(col,'ylabel');
% ylpos=get(yl,'position');
% ylpos(1)=ylpos(1)+2;
% set(yl,'string','Pressure [db]','rotation',270,'position',ylpos)

print('-dpng','-r400',['figures/paper_compare_pressure_',num2str(run1),'_',num2str(run2)])

