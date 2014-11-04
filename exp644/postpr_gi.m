% postprocess gamma^i exp024
clear all
close all

%params; % load some parameters
gamma_i_exp=24;
load(['/home/nfs/z3439823/mymatlab/gamma_i/gamma_i_utils/exp0',num2str(gamma_i_exp),'/data/input_data.mat'])
%keyboard
lat=lats;
lon=longs;

[nz,ny,nx] = size(s);
%keyboard
ibb=backbone_index(squeeze(lon(1,:,:)),squeeze(lat(1,:,:)));

% loading 'values2','vdiff2','gamma_i'
load(['/home/nfs/z3439823/mymatlab/gamma_i/gamma_i_utils/exp0',num2str(gamma_i_exp),'/data/plots_error_3d.mat'])

gi=gamma_i(:,ibb);
pbbi=p(:,ibb);
pbb=nan*values2;
for ii=1:length(values2)
    pbb(ii)=var_on_surf_stef(pbbi,gi,values2(ii));
end
pbb=pbb(~isnan(pbb));
values2_omega=values2(~isnan(pbb));


pbb=pbb(85); %

%keyboard
point=[pbb ibb];
display('optimizing density surface');
%[sns,ctns,pns,rmsdrho,mdf,df_med] =optimize_surface_at_point(s,ct,p,point);
[sns,ctns,pns] =optimize_surface_at_point(s,ct,p,lon,lat,point);
% keyboard
% [sns3d,ctns3d,pns3d,rmsdrho,mdf,df_med] = ribs(ilat,ilon,pbb,s,ct,p);
% 
% save('data/gamma_i_ribs.mat','pns3d','rmsdrho','mdf','df_med','values2_omega')
% save_netcdf(pns3d,'pns3d','data/gamma_i_pns3d.nc')










