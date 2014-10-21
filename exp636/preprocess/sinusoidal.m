clear all
close all
%restoredefaultpath
%addpath(genpath('../../../../gsw_matlab_v3_02'))

nx=80;
ny=40;
nz=101;

x=linspace(0,2*pi*4,nx);
y=linspace(0,2*pi*4/floor(nx/ny),ny);
[X,Y]=meshgrid(x,y);
Z1=sin(X).*sin(Y);

x=linspace(0,2*pi*2,nx);
y=linspace(0,2*pi*2/floor(nx/ny),ny);
[X,Y]=meshgrid(x,y);
Z2=sin(X).*sin(Y);

s=linspace(35,37,nz);
t=linspace(15,5,nz);
p=linspace(0,4000,nz);

t=repmat(t, [nx 1]);
t=repmat(permute(t,[3 2 1]),[ny 1 1]);
t=permute(t,[2 1 3]);
s=repmat(s, [nx 1]);
s=repmat(permute(s,[3 2 1]),[ny 1 1]);
s=permute(s,[2 1 3]);
p=repmat(p, [nx 1]);
p=repmat(permute(p,[3 2 1]),[ny 1 1]);
p=permute(p,[2 1 3]);

Z1=repmat(permute(Z1,[3 1 2]),[nz 1 1]);
Z2=repmat(permute(Z2,[3 1 2]),[nz 1 1]);

%%%%%%%%%%%%%%%%%%%%
% exponential function
[XX,YY]=meshgrid(0:nx-1,0:ny-1);
cx=floor(nx/2);
cy=floor(ny/2);
ex=exp(-(1/(10*nx))*((XX-cx).^2+(nx/ny)*(YY-cy).^2));

contourf(ex)
colorbar
ex=repmat(permute(ex,[3 1 2]),[nz 1 1]);

%%%%%%%%%%%%%%%%%%%%

t=t.*(1+0.2*Z1)+20*ex;
s=s.*(1+0.02*Z2);

r=gsw_rho(s,t,0*p);

lat=linspace(-80,80,ny);
lon=linspace(0,360,nx+1);
lon=lon(1:end-1);

lat=repmat(lat, [nx 1]);
lats=repmat(permute(lat,[3 2 1]),[nz 1 1]);

lon=repmat(lon, [ny 1]);
longs=repmat(permute(lon,[3 1 2]),[nz 1 1]);

figure
contourf(squeeze(r(1,:,:)))
colorbar
figure
contourf(squeeze(r(:,8,:)))
colorbar


ct=t;
vars = {'s','ct','p','lats','longs'};
save('data/input_data.mat',vars{:})
