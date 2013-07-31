
s1=36;
ct1=10.;
p1=0;

s2=s1;
ct2=ct1+5;
p2=2000;

SA=[s1 s2];

CT=[ct1 ct2];
p=[p1 p2];
lat=[45. 45.];

[N2, N2_p, N2_rho, N2_alpha ,N2_beta]=gsw_Nsquared_min(SA,CT,p,lat);
SA_out = gsw_stabilise_SA_neutral(SA,CT,p);
s1=SA_out(1);
s2=SA_out(2);

s3=34.16;
s3=34.44;
s3=33.8565;
ct3=5;
ct3=6.56;
ct3=3.08;
p3=0.5*(p1+p2);

%[SAns,CTns,pns] = depth_ntp(s3,ct3,p3,[s1 s2]',[ct1 ct2]',[p1 p2]')
[SAns,CTns,pns] = depth_ntp_orig(s3,ct3,p3,[s1 s2]',[ct1 ct2]',[p1 p2]')
