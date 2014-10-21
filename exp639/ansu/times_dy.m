function f=times_dy(f)

load('data/dxdy.mat')

dy=0.5*(dy+circshift(dy,[1 0]));
dy(1,:)=dy(2,:); % sloppy

f=f.*dy;

