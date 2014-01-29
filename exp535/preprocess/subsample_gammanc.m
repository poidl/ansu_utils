clear all

lat=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','lat');
lon=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','lon');
p=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','pressure');
s=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','s');
t=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','t');
gamma=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','gamma');

lat=lat(2:end-1); % cut what is not in gk_ak_gamma.mat
s=s(:,:,2:end-1);
t=t(:,:,2:end-1);
gamma=gamma(:,:,2:end-1);

lat=repmat(lat,[1 90]); lat=repmat(lat,[1 1 33 ]);
lon=repmat(lon,[1 43]); lon=repmat(lon,[1 1 33 ]);
lat=permute(lat,[3 1 2]);
lon=permute(lon,[3 2 1]);

p=repmat(p,[1,90]); p=repmat(p,[1 1 43]);

s=permute(s, [1 3 2]);
t=permute(t, [1 3 2]);
p=permute(p, [1 3 2]);

% convert sp to sa
s=gsw_SA_from_SP(s,p,lon,lat);

ct=gsw_CT_from_t(s,t,p);

lats=lat;
longs=lon;

s(s==0)=nan;

vars = {'s','ct','p','lats','longs'};
save('data/input_data.mat',vars{:})




