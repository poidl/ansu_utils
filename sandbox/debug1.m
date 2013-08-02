restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../../../stabilization_paul'))

close all;
clear all;


s3=33.4552;
ct3=-1.7425;
p3=29.1309;

load cast.mat
s12=sc;
ct12=ctc;
p12=pc;

ii=2;
s12=s12(1:ii);
ct12=ct12(1:ii);
p12=p12(1:ii);

s3_=s3*ones(size(p12));
ct3_=ct3*ones(size(p12));
pref=0.5*(p12+p3);
F=gsw_rho(s12,ct12,pref)-gsw_rho(s3_,ct3_,pref);

plot(F,-p12,'*')
hold on
plot([0 0],get(gca,'ylim')) 
hold off

disp('Method             SA       CT       P')
[SAns,CTns,pns] = depth_ntp_jackett(s3,ct3,p3,s12',ct12',p12');
disp(['Jackett:           ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])
[SAns,CTns,pns] = depth_ntp_guillaume(s3,ct3,p3,s12',ct12',p12');
disp(['Guillaume:         ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])
[SAns,CTns,pns] = depth_ntp_jackett_fzero(s3,ct3,p3,s12',ct12',p12');
disp(['Jackett   (fzero): ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])
%[SAns,CTns,pns] = depth_ntp_guillaume_fzero(s3,ct3,p3,s12',ct12',p12');
%disp(['Guillaume (fzero): ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])
[SAns,CTns,pns] = depth_ntp_iter(s3,ct3,p3,s12,ct12,p12);
disp(['depth_ntp_iter:    ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])

