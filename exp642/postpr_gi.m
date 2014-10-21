% postprocess gamma^i exp024
clear all
close all

%params; % load some parameters
gamma_i_exp=24;
load(['/home/nfs/z3439823/mymatlab/gamma_i/gamma_i_utils/exp0',num2str(gamma_i_exp),'/data/input_data.mat'])

% the lats/longs are only used to calculate epsilon in diagnose_and_write()
% they are not necessary to calculate omega surfaces
lat=squeeze(lats(1,:,:)); lon=squeeze(longs(1,:,:));
[dy,dx]=scale_fac(lat,lon);
save('data/dxdy.mat', 'dx', 'dy') 
save('data/latlon.mat', 'lat', 'lon')

[nz,ny,nx] = size(s);

lat=lats;
lon=longs;
% boundary
bdy= 187<=lon(:) & lon(:)<=189 & -17<=lat(:) & lat(:)<=-15; % 16 S, 188 E
if sum(bdy)~=nz
    disp('WARNING: multiple backbones (or none at all)')
    keyboard
end
lon1=squeeze(lon(1,1,:));
lat1=squeeze(lat(1,:,1));
blon=187<=lon1 & lon1<=189; % 16 S, 188 E
blat=-17<=lat1 & lat1<=-15;

if (sum(blon)~=1 || sum(blat)~=1)
    error('something is wrong')
end
ilon=find(blon);
ilat=find(blat);

istation=ilat+ny*(ilon-1);
save('data/stationindex.mat','istation')

% loading 'values2','vdiff2','gamma_i'
load(['/home/nfs/z3439823/mymatlab/gamma_i/gamma_i_utils/exp0',num2str(gamma_i_exp),'/data/plots_error_3d.mat'])

gi=gamma_i(:,ilat,ilon);
pbbi=p(:,ilat,ilon);
pbb=nan*values2;
for ii=1:length(values2)
    pbb(ii)=var_on_surf_stef(pbbi,gi,values2(ii));
end
pbb=pbb(~isnan(pbb));
values2_omega=values2(~isnan(pbb));

pbb=pbb(85); %

%keyboard
point=[pbb ilat ilon];
display('optimizing density surface');
%[sns,ctns,pns,rmsdrho,mdf,df_med] =optimize_surface_at_point(s,ct,p,point);
[sns,ctns,pns] =optimize_surface_at_point(s,ct,p,point);
% keyboard
% [sns3d,ctns3d,pns3d,rmsdrho,mdf,df_med] = ribs(ilat,ilon,pbb,s,ct,p);
% 
% save('data/gamma_i_ribs.mat','pns3d','rmsdrho','mdf','df_med','values2_omega')
% save_netcdf(pns3d,'pns3d','data/gamma_i_pns3d.nc')










