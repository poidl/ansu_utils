clear all
close all

%params; % load some parameters

load('data/input_data.mat')

% clean up
kinds=1:size(s,1); 
for ii=1:size(s,3)
    for jj=1:size(s,2)
        igood=~isnan(s(:,jj,ii));
        if any(igood)
            kg=find(igood,1,'first'); % kg may not be equal to zero: sea ice
            deeper=(kinds' > kg);
            kk=find(isnan(s(:,jj,ii)) & deeper,1,'first'); % shallowest nan below data
            s(kk:end,jj,ii)=nan;
            ct(kk:end,jj,ii)=nan;
        end
    end
end

% the lats/longs are only used to calculate epsilon in diagnose_and_write()
% they are not necessary to calculate omega surfaces
lat=squeeze(lats(1,:,:)); lon=squeeze(longs(1,:,:));
[dy,dx]=scale_fac(lat,lon);
save('data/dxdy.mat', 'dx', 'dy') 
save('data/latlon.mat', 'lat', 'lon') 

sa=s; clear s; % the _subs_ data saves sa in variable 's'

[zi,yi,xi] = size(sa);

make_land_mask

%omega_3d

% ilat=18;
% ilon=48;
% dep=1960;

ilat=5;
ilon=85;
dep=36.748887365683913;

% ilat=29;
% ilon=69;
% dep=532.2648;

istation=ilat+yi*(ilon-1);
save('data/stationindex.mat','istation')

load('../exp603/data/iteration_history.mat')
sns=squeeze(sns_hist(end,:,:));
ctns=squeeze(ctns_hist(end,:,:));
pns=squeeze(pns_hist(end,:,:));

[sns,ctns,pns] = optimize_surface_exact(sa,ct,p,sns,ctns,pns);




        
        
  