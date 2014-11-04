function [dx,dy,dz]=get_dx(lon,lat,p)

la=squeeze(lat(1,:,:));
lo=squeeze(lon(1,:,:));

[dy,dx]=scale_fac(la,lo);
%save('data/dxdy.mat', 'dx','dy') 
%load('data/dxdy.mat');
[nz,ny,nx]=size(lat);
dx=repmat(permute(dx,[3 1 2]),[nz 1 1]);
dy=repmat(permute(dy,[3 1 2]),[nz 1 1]);
dz=circshift(p,[-1 0 0])-p;
dz(end,:,:)=nan;

% dx=0*dx+1;
% dy=0*dy+1;
% dz=0*dz+1;