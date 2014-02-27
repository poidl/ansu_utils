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

sz=0.8*[25 21];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)])


subplot(2,2,1)
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


subplot(2,2,3)
h2=imagesc(lon,lat,vv1it);
set(gca,'YDir','normal')
set(h2,'alphadata',~isnan(vv1it))
cb2=colorbar()
hold on
plot(lon(iis),lat(jis),'kx','markersize',20,'linewidth',3)

xlabel('Longitude')
ylabel('Latitude')
load ('../external_scripts/coast/coast_data.mat');
plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot b on the surface

subplot(2,2,2)
run1=571;
fname=['../exp',num2str(run1),'/data/iteration_history.mat'];
varname= 'b_hist';
load(fname, varname);

vv1=b_hist; % variable to plot

vv2it=squeeze(vv1(it-1,:,:));

h3=imagesc(lon,lat,vv2it);
set(gca,'YDir','normal')
set(h3,'alphadata',~isnan(vv2it))
cb3=colorbar()
hold on
plot(lon(iis),lat(jis),'kx','markersize',20,'linewidth',3)

xlabel('Longitude')
ylabel('Latitude')
load ('../external_scripts/coast/coast_data.mat');
plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);


disp(['min b:',num2str(min(vv2it(:)))])
disp(['max b:',num2str(max(vv2it(:)))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmp=colormap(hot(128));
cmp2=cmp(find(cmp(:,1)==1,1,'first')-15:end,:);
cm1=colormap([fliplr(cmp2);flipud(cmp2)]) ;
cm2=jet;
cmap=[cm1;cm2;cm1];

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

va=vv2it;
cbb=max(abs(va(:)-1));
cbb=0.25;
cmax3=1+cbb;
cmin3=1-cbb;
C3 = min(m,round((m-1)*(va-cmin3)/(cmax3-cmin3))+1);
C3 = 2*m+C3;

set(h1,'CData',C1);
set(h2,'CData',C2);
set(h3,'CData',C3);

ax = findobj(gcf,'Type','axes');
set(ax,'CLim', [1 3*m])

set(cb1,'ylim',[1 m])
set(cb2,'ylim',[m+1 2*m])
set(cb3,'ylim',[2*m+1 3*m])

nticks=5;
set(cb1,'ytick',linspace(1,m,nticks))
set(cb1,'yticklabel',linspace(cmin1,cmax1,nticks))

set(cb2,'ytick',linspace(m+1,2*m,nticks))
set(cb2,'yticklabel',linspace(cmin,cmax,nticks))

set(cb3,'ytick',linspace(2*m+1,3*m,nticks))
set(cb3,'yticklabel',linspace(cmin3,cmax3,nticks))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% position
width=1.4;
left=0.06;
cbwidth=0.7;
cbshift=-0.005;
lshift=0.6;
%%%
p1=get(ax(6),'position');
p1(3)=width*p1(3);
p1(1)=left;
yp=get(cb1,'position');
yp(3)=cbwidth*yp(3);
yp(1)=cbshift+yp(1);
set(cb1,'position',yp);
set(ax(6),'position',p1);

yl=get(cb1,'ylabel');
ylpos=get(yl,'position');
ylpos(1)=ylpos(1)+lshift;
%set(yl,'position',ylpos)
set(yl,'string','Pressure difference EXP1''-EXP5 [db]','rotation',270,'position',ylpos)
% 
% %%%
p1=get(ax(4),'position');
p1(3)=width*p1(3);
p1(1)=left;
yp=get(cb2,'position');
yp(3)=cbwidth*yp(3);
yp(1)=cbshift+yp(1);
set(cb2,'position',yp);
set(ax(4),'position',p1);

yl=get(cb2,'ylabel');
ylpos=get(yl,'position');
ylpos(1)=ylpos(1)+lshift;
%set(yl,'position',ylpos)
set(yl,'string','Pressure EXP1'' [db]','rotation',270,'position',ylpos)
% %%%
p1=get(ax(2),'position');
p1(3)=width*p1(3);
p1(1)=left+0.52;
yp=get(cb3,'position');
yp(3)=cbwidth*yp(3);
yp(1)=cbshift+yp(1)+0.08;
set(cb3,'position',yp);
set(ax(2),'position',p1);
 
yl=get(cb3,'ylabel');
ylpos=get(yl,'position');
ylpos(1)=ylpos(1)+lshift;
%set(yl,'position',ylpos)
set(yl,'string','b^*','rotation',270,'position',ylpos)

text(-0.15,1.1,'b)','units','normalized','fontsize',15,'Parent',ax(2))
text(-0.15,1.1,'c)','units','normalized','fontsize',15,'Parent',ax(4))
text(-0.15,1.1,'a)','units','normalized','fontsize',15,'Parent',ax(6))


print('-dpng','-r400',['figures/paper_compare_pressure_',num2str(run1),'_',num2str(run2)])

