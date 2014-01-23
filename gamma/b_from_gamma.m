clear all
close all
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))

lat=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','lat');
lon=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','lon');
p=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','pressure');
s=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','s');
t=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','t');
gamma=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','gamma');

p=repmat(p,[1,90]); p=repmat(p,[1 1 45]);

ct=gsw_CT_from_t(s,t,p);

smid=s(1:end-1,:,:)+0.5*diff(s,1,1);
ctmid=ct(1:end-1,:,:)+0.5*diff(ct,1,1);
pmid=p(1:end-1,:,:)+0.5*diff(p,1,1);

rho=gsw_rho_CT(smid,ctmid,pmid);

[n2,pmid]=gsw_Nsquared(s,ct,p);

n2=reshape(n2,[32 90 45]);
pmid=reshape(pmid,[32 90 45]);

n3=(-rho./9.81).*n2;

dgamma=-(gamma(1:end-1,:,:)-gamma(2:end,:,:))./(p(1:end-1,:,:)-p(2:end,:,:));

b=dgamma./n3;

b=permute(b, [1 3 2]);

p=permute(pmid, [1 3 2]);
lat=repmat(lat,[1 90]); lat=repmat(lat,[1 1 32 ]);
lon=repmat(lon,[1 45]); lon=repmat(lon,[1 1 32 ]);
lat=permute(lat,[3 1 2]);
lon=permute(lon,[3 2 1]);


save b_original.mat lat lon p b








