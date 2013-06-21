clear all
close all

user_input; % load configuration from settings.m

load('data/input_data.mat')
lats=squeeze(lats(1,:,1)); longs=squeeze(longs(1,1,:))'; p=squeeze(p(:,1,1));
sa=s; clear s; % the _subs_ data saves sa in variable 's'


%g = gsw_grav(lats);

[zi,yi,xi] = size(sa);


%[n2tmp,p_mid] = gsw_Nsquared(sa,ct,p,lats);
%n2=reshape(n2tmp,size(sa)-[1 0 0]);

rho = gsw_rho(sa,ct,zeros(size(sa)))-1e3;

lon=200;
ilon=find(longs>=lon,1,'first');
vv=squeeze(rho(:,:,ilon))
cf=contour(lats,-p,vv);
colorbar
hold on 
contour(lats,-p,vv, 27.7*[1 1])
%set(gca,'YDir','normal')