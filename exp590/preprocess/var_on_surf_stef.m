function vsurf=var_on_surf_stef(va,p,surf)
% variable is 3d and defined on p
% surf is surface of p

[nz,ny,nx]=size(va);

p0=surf(:)';
p_stacked=repmat(p0, [nz 1]);

up=p(:,:)<=p_stacked;
kup=sum(up,1); % number of grid points above surface

inan=isnan(p0);
bottom= kup==nz;

kup=kup+nz*(0:nx*ny-1); % flat index

kup(inan|bottom)=1; % dummy

va_1=va(kup);
va_2=va(kup+1);
p_1=p(kup);
p_2=p(kup+1);

dp=(p0-p_1)./(p_2-p_1);
vsurf=va_1+(va_2-va_1).*dp;
upper_only= dp==0; % only consider upper bottle where dp==0 (in case lower bottle is nan)
vsurf(upper_only)=va_1(upper_only);

vsurf(inan)=nan;
vsurf(bottom)=va(end,bottom);

vsurf=reshape(vsurf,[ny nx]);

end