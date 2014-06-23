close all
clear all

load('data/bb_coords_surface.mat')
dlat=diff(ilat_);
dlat=dlat~=0;
dlon=diff(ilon_);
dlon=dlon~=0;
dl=dlat|dlon;
ns_surf=sum(dl)+1;

load('data/bb_coords_bottom.mat')
dlat=diff(ilat_);
dlat=dlat~=0;
dlon=diff(ilon_);
dlon=dlon~=0;
dl=dlat|dlon;
ns_bot=sum(dl)+1;

nk_surface=[];
nk_bottom=[];
ps=[];

for ii=ns_bot+1:ns_bot+ns_surf
    
    s1=ncread(['data/nc/surface/sns3d/sns3d',num2str(ii),'.nc'],'sns3d');
    ct1=ncread(['data/nc/surface/ctns3d/ctns3d',num2str(ii),'.nc'],'ctns3d');
    p1=ncread(['data/nc/surface/pns3d/pns3d',num2str(ii),'.nc'],'pns3d');
    s1=permute(s1,[3 2 1]);
    ct1=permute(ct1,[3 2 1]);
    p1=permute(p1,[3 2 1]);    
    if isempty(ps)
        ss=s1;
        cts=ct1;        
        ps=p1;
    else
        ss=cat(1,ss,s1);
        cts=cat(1,cts,ct1);
        ps=cat(1,ps,p1);
    end
    nk_surface=[nk_surface,size(p1,1)];
end

for ii=1:ns_bot
    s1=ncread(['data/nc/bottom/sns3d/sns3d',num2str(ii),'.nc'],'sns3d');
    ct1=ncread(['data/nc/bottom/ctns3d/ctns3d',num2str(ii),'.nc'],'ctns3d');
    p1=ncread(['data/nc/bottom/pns3d/pns3d',num2str(ii),'.nc'],'pns3d');
    s1=permute(s1,[3 2 1]);
    ct1=permute(ct1,[3 2 1]);
    p1=permute(p1,[3 2 1]); 

    ss=cat(1,ss,s1);
    cts=cat(1,cts,ct1);
    ps=cat(1,ps,p1);    
    
    nk_bottom=[nk_bottom,size(p1,1)];
end

nk_surface=cumsum(nk_surface);
nk_bottom=cumsum(nk_bottom);

load('data/bb_coords_surface.mat')
ilon_s=ilon_;
ilat_s=ilat_;
load('data/bb_coords_bottom.mat')
ilon=[ilon_s,ilon_];
ilat=[ilat_s,ilat_];
load('data/p_bb_surface.mat')
p_bb_s=p_bb;
load('data/p_bb_bottom.mat')
pbb=[p_bb_s,p_bb];

[ilat,ilon,sns3d,ctns3d,pns3d]=sortit(ilat,ilon,ss,cts,ps);
disp('final test: ')
if checkit(pns3d)
    error('not good')
else
    disp('success')
end

save('data/bb_coords_final.mat','ilat','ilon')
save_netcdf(sns3d,'sns3d',['data/nc/sns3d_final.nc'])
save_netcdf(ctns3d,'ctns3d',['data/nc/ctns3d_final.nc'])
save_netcdf(pns3d,'pns3d',['data/nc/pns3d_final.nc'])

keyboard
