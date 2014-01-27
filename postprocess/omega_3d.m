clear all;
close all;

fname=['../exp344/data/stationindex.mat'];
load(fname);
disp(istation)

fname=['../exp344/data/input_data.mat'];
load(fname);

grid_s_i=s(:,istation);
grid_ct_i=ct(:,istation);
grid_p_i=p(:,istation);

load skeleton.mat 

pmid=0.5*(grid_p_i(1:end-1)+grid_p_i(2:end));
r1=gsw_rho(grid_s_i(1:end-1),grid_ct_i(1:end-1),pmid);
r2=gsw_rho(grid_s_i(2:end),grid_ct_i(2:end),pmid);

rgN2=r1-r2;

if grid_p_i(1)~=0
    error('shallowest point not on surface')
end
grid_omega_i=gsw_rho(grid_s_i(1),grid_ct_i(1),grid_p_i(1))-cumsum(rgN2);

nsurf=length(p_istation);
om_istation=nan*ones(1,nsurf);

for ii=1:nsurf
    
    surf_p=p_istation(ii);
    surf_s=s_istation(ii);
    surf_ct=ct_istation(ii);
    
    [mini,ig]=min(abs(grid_p_i<surf_p));
    
    pu=grid_p_i(ig);
    su=grid_s_i(ig);
    ctu=grid_ct_i(ig);
    
    t1=grid_omega_i(ig);
        
    pmid=0.5*(pu+surf_p);
    r1=gsw_rho(su,ctu,pmid);
    r2=gsw_rho(surf_s,surf_ct,pmid);
    
    t2=r1-r2;

    if grid_p_i(ig)<surf_p   
        om_istation(ii)=t1+t2;
    else
        om_istation(ii)=t1-t2;
    end
    
end

d=3;


    




% 
% vp=squeeze(surfvar(1,:,:));
% h=imagesc( vp )
% set(h,'alphadata',~isnan(vp)) % white nans
% set(gca,'YDir','normal')
% colorbar()

