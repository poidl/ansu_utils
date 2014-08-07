function vsurf=var_on_surf_stef(va,p,surf)
% variable is 3d and defined on p
% surf is surface of p

[nz,ny,nx]=size(va);
p0=surf(:)';
p_stacked=repmat(p0, [nz 1]);

up=p(:,:)<=p_stacked;

% don't assume p is withouth nans
% don't assume p is monotonous (in case it's not pressure)

tmp= up & circshift(~up,-1); 
tmp(end,:)=false; % discard bottom
cs=cumsum(tmp,1)+1;
cs(cs~=1)=0;
kup=sum(cs,1)+1;
nothing=~any(tmp,1);

kup=kup+nz*(0:nx*ny-1); % flat index

kup(nothing)=1;

va_1=va(kup);
va_2=va(kup+1);
p_1=p(kup);
p_2=p(kup+1);

dp=(p0-p_1)./(p_2-p_1);
vsurf=va_1+(va_2-va_1).*dp;
upper_only= dp==0; % only consider upper bottle where dp==0 (in case lower bottle is nan)
vsurf(upper_only)=va_1(upper_only);

vsurf(nothing)=nan;

bottom_good=p(end,:)==p_stacked(end,:);
vsurf(bottom_good)=va(end,bottom_good); % bottom exact

vsurf=reshape(vsurf,[ny nx]);

end