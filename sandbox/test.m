restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))

s=linspace(0,39,50);
ct=linspace(-5,25,50);
[s,ct]=meshgrid(s,ct);

rho=gsw_rho(s,ct,0);
contour(rho)
title('pref=0')

figure()
rho=gsw_rho(s,ct,1000);
contour(rho)
title('pref=1000')