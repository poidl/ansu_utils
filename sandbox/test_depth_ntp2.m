% the s,ct,p cast generated here is plotted in stabilize5.m
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../../../stabilization_paul'))

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

[SAns,CTns,pns] = depth_ntp(s3,ct3,p3,[s1 s2]',[ct1 ct2]',[p1 p2]')
[SAns,CTns,pns] = depth_ntp_orig(s3,ct3,p3,[s1 s2]',[ct1 ct2]',[p1 p2]')
[SAns,CTns,pns] = depth_ntp_jackett(s3,ct3,p3,[s1 s2]',[ct1 ct2]',[p1 p2]')
[SAns,CTns,pns] = depth_ntp_iter(s3,ct3,p3,[s1 s2]',[ct1 ct2]',[p1 p2]')
