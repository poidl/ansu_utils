clear all;
close all;

fname=['../exp344/data/stationindex.mat'];
load(fname);
disp(istation)

fname=['../exp344/data/input_data.mat'];
load(fname);

grid_si=s(:,istation);
grid_cti=ct(:,istation);
grid_pi=p(:,istation);

load skeleton.mat 

surf_pi=p_istation;
surf_si=s_istation;
surf_cti=ct_istation;

midgrid_pi=0.5*(grid_pi(1:end-1)+grid_pi(2:end));
r1=gsw_rho(grid_si(1:end-1),grid_cti(1:end-1),midgrid_pi);
r2=gsw_rho(grid_si(2:end),grid_cti(2:end),midgrid_pi);

midgrid_rgN2=r1-r2;

if grid_pi(1)~=0
    error('shallowest point not on surface')
end
midgrid_omegai=gsw_rho(grid_si(1),grid_cti(1),grid_pi(1))-cumsum(midgrid_rgN2);

nsurf=length(surf_pi);
surf_omegai=nan*ones(1,nsurf);

midgrid_si=0.5*(grid_si(1:end-1)+grid_si(2:end));
midgrid_cti=0.5*(grid_cti(1:end-1)+grid_cti(2:end));

for ii=1:nsurf
    
    surf_p=surf_pi(ii);
    surf_s=surf_si(ii);
    surf_ct=surf_cti(ii);
    
    [mini,img]=min(abs(midgrid_pi-surf_p));
    img=find(midgrid_pi<surf_p,1,'last');
    
    pu=midgrid_pi(img);
    su=midgrid_si(img);
    ctu=midgrid_cti(img);
    
    t1=midgrid_omegai(img);
        
    pm=0.5*(pu+surf_p);
    r1=gsw_rho(su,ctu,pm);
    r2=gsw_rho(surf_s,surf_ct,pm);
    
    t2=r1-r2; 
    
    surf_omegai(ii)=t1-t2;
    
end



[zi,yi,xi]=size(s);
s=s(:,:);
ct=ct(:,:);
p=p(:,:);

grid_omega=nan*ones(zi,yi*xi);

for ii=1:yi*xi %include ii==istation because so far we only have omega on midgrid at istation
    surf_p=pns_3d(:,ii);
    if ~isnan(sns_3d(10,ii))
        dummy=10;
    end
    for kk=1:zi
        above=find(surf_p<=p(kk,ii),1,'last');
        below=find(surf_p>p(kk,ii),1,'first');
        if ~isempty(above) && ~isempty(below) && below==(above+1)
            dl=(p(kk,ii)-surf_p(above)) / (surf_p(below)-surf_p(above));
            
            grid_omega(kk,ii)=surf_omegai(above)+dl*(surf_omegai(below)-surf_omegai(above));
            
%             pp_station=surf_pi(above)+dl*(surf_pi(below)-surf_pi(above));
% 
%             %[mini,ig]=min(abs(p(:,ii)-pp_station));
%             ig=find(p(:,ii)<pp_station,1,'last');
%             dlg=(pp_station-p(ig,ii))/(p(ig+1,ii)-p(ig,ii));
% 
%             ss_station=grid_si(ig,ii)+dlg*(grid_si(ig+1,ii)-grid_si(ig,ii));
%             ct_station=grid_cti(ig,ii)+dlg*(grid_si(ig+1,ii)-grid_si(ig,ii));
%      
%             %[mini,is]=min(abs(surf_pi-pp_station));
%             is=find(surf_pi<pp_station,1,'last');
%             
%             pu=surf_pi(is);
%             su=surf_si(is);
%             ctu=surf_cti(is);
% 
%             t1=surf_omegai(is);
% 
%             pm=0.5*(pu+pp_station);
%             r1=gsw_rho(su,ctu,pm);
%             r2=gsw_rho(ss_station,ct_station,pm);
% 
%             t2=r1-r2;
%   
%             grid_omega(kk,ii)=t1-t2;
        end
    end
end



omega=reshape(grid_omega,[zi,yi,xi]);

save omega.mat omega

delete omega.nc        
nccreate('omega.nc','omega',...
         'Dimensions', {'k', zi, 'j', yi, 'i', xi})
ncwrite('omega.nc','omega', omega);

% 
% % 
% % vp=squeeze(surfvar(1,:,:));
% % h=imagesc( vp )
% % set(h,'alphadata',~isnan(vp)) % white nans
% % set(gca,'YDir','normal')
% % colorbar()
% 
