function [ex,ey]=times_sqrtdA_on_delta(ex,ey,dx,dy)

dxt=regrid_new(dx,2,-2); % dx onto dy
dx_=regrid_new(dxt,1,2);
dyt=regrid_new(dy,1,-2); % dy onto dx
dy_=regrid_new(dyt,2,2);
f1=sqrt(dx.*dy_);
f2=sqrt(dx_.*dy); 
ex=f1.*ex;
ey=f2.*ey;
