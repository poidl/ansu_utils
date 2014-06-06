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


omega_3d

% lat=lat(:,1);
% lon=lon(1,:);
% [mini,ilat]=min(abs(lat+16)); % Jackett & McDougall 97: 16 South 188 East
% [mini,ilon]=min(abs(lon-188));
% [zi,yi,xi]=size(sa);
% istation=ilat+yi*(ilon-1);
% save('data/stationindex.mat','istation')


%save_netcdf(sa,ct,p,sns,ctns,pns);

% point=[1000 ilat ilon];
% 
% display('optimizing density surface');
% tic
% [sns,ctns,pns] =optimize_surface_at_point(sa,ct,p,point);
% display(['optimizing density surface took ',num2str(toc),' seconds']);
% 
% save('data/surface_final.mat','sns','ctns','pns')

% figure()
% h=imagesc(pns);
% set(h,'alphadata',~isnan(pns));
% set(gca,'YDir','normal');


        
        
  
