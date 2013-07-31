restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../../../stabilization_paul'))

close all;
clear all;

s1=37;
ct1=14.;
p1=0;

s2=36;
ct2=13.9;
p2=2000;

SA=[s1 s2];

CT=[ct1 ct2];
p=[p1 p2];
lat=[45. 45.];

[N2, N2_p, N2_rho, N2_alpha ,N2_beta]=gsw_Nsquared_min(SA,CT,p,lat);
SA_out = gsw_stabilise_SA_neutral(SA,CT,p);
s1=SA_out(1);
s2=SA_out(2);

s3=36.4274;
ct3=13.709;
p3=0.5*(p1+p2);

s12=linspace(s1,s2,100);
ct12=linspace(ct1,ct2,100);
p12=linspace(p1,p2,100);

p13=0.5*(p1+p3);
p23=0.5*(p2+p3);

s3_=s3*ones(size(p12));
ct3_=ct3*ones(size(p12));
pref=0.5*(p12+p3);
F=gsw_rho(s12,ct12,pref)-gsw_rho(s3_,ct3_,pref);

plot(F,-p12)
hold on
plot([0 0],get(gca,'ylim')) 
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_=linspace(36.4,36.55,50);
ct_=linspace(13.65,14.1,50);
[s_,ct_]=meshgrid(s_,ct_);

rho1=gsw_rho(s_,ct_,p1);
rho2=gsw_rho(s_,ct_,p2);
rho13=gsw_rho(s_,ct_,p13);
rho23=gsw_rho(s_,ct_,p23);


figure()
%contour(s_,ct_,rho1)
%hold on
contour(s_,ct_,rho1,gsw_rho([s1 s1],[ct1 ct1],[p1 p1]),'color','k');
hold on
contour(s_,ct_,rho2,gsw_rho([s2 s2],[ct2 ct2],[p2 p2]),'color','k');
contour(s_,ct_,rho13,gsw_rho([s1 s1],[ct1 ct1],[p13 p13]),'color','k','linestyle','--');
contour(s_,ct_,rho23,gsw_rho([s2 s2],[ct2 ct2],[p23 p23]),'color','k','linestyle','--');
plot(s1,ct1,'r*')
plot(s2,ct2,'b*')
plot(s3,ct3,'k*')
plot(s12,ct12)
hold off
title('pref=0')



% figure()
% rho=gsw_rho(s_,ct_,1000);
% contour(rho)
% title('pref=1000')