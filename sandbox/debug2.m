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


s3=33.4552;
ct3=-1.7425;
p3=29.1309;

load cast.mat
s12=[33 34];
ct12=[-1 -2];
p12=[0 100];


s3_=s3*ones(size(p12));
ct3_=ct3*ones(size(p12));
pref=0.5*(p12+p3);
F=gsw_rho(s12,ct12,pref)-gsw_rho(s3_,ct3_,pref);

plot(F,-p12)
hold on
plot([0 0],get(gca,'ylim')) 
hold off

[a,b,c]=depth_ntp_orig(s3,ct3,p3,[33 34]',[-1 -2]',[0 100]')
[a,b,c]=depth_ntp(s3,ct3,p3,[33 34]',[-1 -2]',[0 100]')

[a,b,c]=depth_ntp_orig(s3,ct3,p3,s12',ct12',p12')
[a,b,c]=depth_ntp(s3,ct3,p3,s12,ct12,p12)
