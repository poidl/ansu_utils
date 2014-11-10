clear all
close all
restoredefaultpath
addpath(genpath('../../../neutral_common'))
addpath(genpath('.'))

path_to_gammanc='/home/nfs/z3439823/mymatlab/neutral_common/gamma_n/gamma_JackettMcDougall96';
[s,ct,p,lon,lat]=get_input_jackett96(path_to_gammanc);
%fname='/home/nfs/z3439823/mymatlab/neutral_common/data/WOA13/woa13.nc';
%[s,ct,p,lon,lat]=get_input_woa13_4deg(fname);
%fname='/home/nfs/z3439823/mymatlab/neutral_common/data/wghc/convert_to_netcdf.py/wghc.nc';
%[s,ct,p,lon,lat]=get_input_wghc_4deg(fname);  
%fname='/home/nfs/z3439823/mymatlab/neutral_common/data/nemo/nemo.nc';
%[s,ct,p,lon,lat]=get_input_nemo_4deg(fname); 
        
[s,ct,p,lon,lat]=crop_to_gamma_n(s,ct,p,lon,lat);
[s,ct]=only_keep_largest_region(s,ct); % remove everything except largest region

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pbb is the vector of pressures at which the surfaces are clamped at the 
% backbone (bb). ibb is the flat lateral index of the backbone location. 
% ibb and the elements of pbb 
% form the 'point' argument for optimize_surface_at_point().
% Below we construct a pbb, but this is not necessary. You can
% find a single surface by running
% optimize_surface_at_point(s,ct,p,lon,lat,[p0 ibb]); 

lon_=squeeze(lon(1,:,:));
lat_=squeeze(lat(1,:,:));
[ibb,ilat,ilon]=backbone_index(lon_,lat_); % open to change lateral location of bb

pbot=p(sum(~isnan(s(:,ibb))),ibb);
pbb=linspace(0,pbot,102);
pbb=pbb(2:end-1);
% get gamma_n values at clamping points
sbb=interp1(p(:,ibb),s(:,ibb),pbb);  
ctbb=interp1(p(:,ibb),ct(:,ibb),pbb);
gvals=get_gamma_n(sbb',ctbb',pbb',lon(:,ibb),lat(:,ibb));

ns=length(pbb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,yn,xn]=size(s);
sns3d=nan*ones(ns,yn,xn);
ctns3d=sns3d;
pns3d=sns3d;
sms=nan*ones(1,ns);
dms=nan*ones(1,ns);

%pbb=pbb(50);
%ns=length(pbb);

for kk=1:ns
    point=[pbb(kk) ibb];
    display(['pressure ',num2str(pbb(kk))]);
    [sns,ctns,pns] =optimize_surface_at_point(s,ct,p,lon,lat,point);
    sms(kk)=error_surface(sns,ctns,pns,s,ct,p,lon,lat,'slope_difference');
    dms(kk)=error_surface(sns,ctns,pns,s,ct,p,lon,lat,'drho_local');
    sns3d(kk,:,:)=sns;
    ctns3d(kk,:,:)=ctns;
    pns3d(kk,:,:)=pns;
end

omega_user_input;
save([datapath,'/omega_3d.mat'],'pns3d','sms','dms')
save_netcdf03(pns3d,'pns3d',[datapath,'/omega_3d.nc'])
