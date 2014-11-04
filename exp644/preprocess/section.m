clear all
close all

omega_user_input; % load configuration from settings.m

load('data/input_data.mat')
lats=squeeze(lats(1,:,1)); longs=squeeze(longs(1,1,:))'; p=squeeze(p(:,1,1));
sa=s; clear s; % the _subs_ data saves sa in variable 's'


%g = gsw_grav(lats);

[zi,yi,xi] = size(sa);


%[n2tmp,p_mid] = gsw_Nsquared(sa,ct,p,lats);
%n2=reshape(n2tmp,size(sa)-[1 0 0]);
mypr=1000;
rho = gsw_rho(sa,ct,mypr*ones(size(sa)))-1e3;
rho_zm=nanmean(rho,3);

if 0
    lon=200;
    ilon=find(longs>=lon,1,'first');
    vv=squeeze(rho(:,:,ilon))
end
vv=squeeze(rho_zm)
%cf=contour(lats,-p,vv);
%hold on
l1=linspace(33,37, 9); %pr 2000
l2=linspace(33,37, 5);
l1=linspace(24,28, 9); %pr 0
l2=linspace(24,28, 5);
l1=linspace(29,33,9); %pr 1000
l2=linspace(29,33, 5);
%l1=linspace(25,29,9); %pr 300
%l2=linspace(25,29, 5);
c=contour(lats,-p,vv,l1,'color','k')
hold on
c=contour(lats,-p,vv,l2,'color','k','linewidth',1)
clabel(c)
% c=contour(lats,-p,vv,'color','k','linewidth',1)
% clabel(c)
title(['Zonal mean of \sigma_{pr=',num2str(mypr),'}'],'fontsize',15)

print('-dpng','-r200',['figures/rho_zm_contours_for_pr_',num2str(mypr)])