close all
clear all

nflat=1;
nslope=1;
ny=2*nflat;
nz=101;
nx=2;

fac=[ones(1,nflat), linspace(1,0,nslope), zeros(1,nflat)];
fac=[ones(1,nflat), zeros(1,nflat)];

theta1=1.8;
theta2=-0.2;
cast1=linspace(theta1,theta2,nz)';
dt=0.4;
cast2=linspace(theta1+dt,theta2+dt,nz)';

t1=bsxfun(@times, ones(1,2*nflat) , cast2);

t2=bsxfun(@times, fac, cast1-cast2);

ct=t1+t2;

imagesc(ct)
colorbar

ct=repmat(ct,[1,1,nx]);

p=linspace(0,5000,nz)';
p=repmat(p,[1,ny]);
p=repmat(p,[1,1,nx]);

lat=linspace(-85,85,ny);
lat=bsxfun(@times,lat,ones(nz,1));
lats=repmat(lat,[1,1,nx]);
lon=linspace(0,4,nx);
lon=bsxfun(@times,lon,ones(ny,1));
longs=repmat(permute(lon,[3 1 2]),[nz 1 1]);

s=35+0*ct;

rho0=gsw_rho(squeeze(s(:,:,1)), squeeze(ct(:,:,1)), 0*squeeze(p(:,:,1)));
if 1;
    figure
    contourf(lats(1,:,1),-p(:,1,1),rho0)
    colorbar
end

vars = {'s','ct','p','lats','longs'};
save('data/idealized_01.mat',vars{:})

