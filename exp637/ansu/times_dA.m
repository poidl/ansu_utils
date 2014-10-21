function f=times_dA(f)

load('data/dxdy.mat')

dx=0.5*(dx+circshift(dx,[0 1]));
dx(:,1)=dx(:,2); % sloppy
dy=0.5*(dy+circshift(dy,[1 0]));
dy(1,:)=dy(2,:); % sloppy

dA=dx.*dy;

f=f.*dA;

