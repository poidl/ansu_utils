function f=times_dx(f)

load('data/dxdy.mat')

dx=0.5*(dx+circshift(dx,[0 1]));
dx(:,1)=dx(:,2); % sloppy

f=f.*dx;
