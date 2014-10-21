function f=over_dA(f,staggering)
% staggering="i": f lives on the same j-grid as s,t,p, but is staggered in
% the i-dimension
% staggering="j": f lives on the same i-grid as s,t,p, but is staggered in
% the j-dimension
load('data/dxdy.mat')

if strcmp(staggering,'i')
    
    dy=0.5*(dy+circshift(dy,[1 0]));
    dy(1,:)=dy(2,:); % sloppy
    dy=0.5*(dy+circshift(dy,[0 -1]));
    dy(:,end)=dy(:,end-1); % sloppy

elseif strcmp(staggering,'j')
    
    dx=0.5*(dy+circshift(dx,[0 1]));
    dx(:,1)=dx(:,2); % sloppy
    dx=0.5*(dx+circshift(dx,[-1 0]));
    dx(end,:)=dx(end-1,:); % sloppy
    
else
    error('problem')
end

dA=dx.*dy;
keyboard
f=f./dA;
