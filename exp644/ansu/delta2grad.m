function [fx,fy]=delta2grad(fx,fy)

load('data/dxdy.mat')

nxy=length(dx(:));

if length(fx(:))==nxy % fx is probably a lateral surface
    fx(:)=fx(:)./dx(:);
    fy(:)=fy(:)./dy(:);
elseif size(fx(:,:),2)==nxy % fx is probably a 3d field
%    keyboard
    fx(:,:)=bsxfun(@times,fx(:,:),dx(:)');
    fy(:,:)=bsxfun(@times,fy(:,:),dy(:)');
else
    error('problem')
end


