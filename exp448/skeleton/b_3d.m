clear all;
close all;


addpath(genpath('../../../../gsw_matlab_v3_02'))

fname=['../../exp344/data/stationindex.mat'];
load(fname);
disp(istation)

fname=['../../exp344/data/input_data.mat'];
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

b=nan*ones(zi,yi*xi);

domega=surf_omegai(1:end-1)-surf_omegai(2:end);
midgrid_domega=repmat(domega',[1,xi*yi]);

sns_3d=sns_3d;
ctns_3d=ctns_3d;
pmid=0.5*(pns_3d(1:end-1,:)+pns_3d(2:end,:));
r1=gsw_rho(sns_3d(1:end-1,:),ctns_3d(1:end-1,:),pmid);
r2=gsw_rho(sns_3d(2:end,:),ctns_3d(2:end,:),pmid);
midgrid_n3=r1-r2;

midgrid_b=midgrid_domega./midgrid_n3;


for ii=1:yi*xi %include ii==istation because so far we only have omega on midgrid at istation
    surf_p=pns_3d(:,ii);
    if ~isnan(sns_3d(10,ii))
        dummy=10;
    end
    for kk=1:zi
        above=find(surf_p<=p(kk,ii),1,'last');
        below=find(surf_p>p(kk,ii),1,'first');
        if ~isempty(above) && ~isempty(below) && below==(above+1)
            
%             domega=surf_omegai(above)-surf_omegai(below);
%             domega=domega*ones(1,yi*xi);
%             
%             pmid=0.5*(pns_3d(above,:)+pns_3d(below,:));
%             r1=gsw_rho(sns_3d(above,:),ctns_3d(above,:),pmid);
%             r2=gsw_rho(sns_3d(below,:),ctns_3d(below,:),pmid);
%             n3=r1-r2;
            
%            b(kk,ii)=domega(ii)./n3(ii);

            b(kk,ii)=midgrid_b(above,ii);
        end
    end
end

% clear all;
% close all;
% 
% fname=['../exp344/data/input_data.mat'];
% load(fname);
% 
% [zi,yi,xi]=size(s);
% 
% p=p(:,:);
% s=s(:,:);
% ct=ct(:,:);
% 
% load omega.mat 
% 
% omega=omega(:,:);
% 
% pmid=0.5*(p(1:end-1,:)+p(2:end,:));
% r1=gsw_rho(s(1:end-1,:),ct(1:end-1,:),pmid);
% r2=gsw_rho(s(2:end,:),ct(2:end,:),pmid);
% 
% dz=-(p(1:end-1,:)-p(2:end,:));
% 
% mid_rgN2=(r1-r2)./dz;
% 
% om_z=(omega(1:end-1,:)-omega(2:end,:))./dz;
% 
% b=om_z./mid_rgN2;
% 
b=reshape(b,[zi,yi,xi]);
% pb=reshape(pmid,[zi-1,yi,xi]);
% 
save('../data/b.mat','b') 

delete b.nc        
nccreate('b.nc','b',...
         'Dimensions', {'k', zi, 'j', yi, 'i', xi})
ncwrite('b.nc','b', b);
