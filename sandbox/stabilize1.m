restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../../../stabilization_paul'))

close all;
clear all;

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

s12=linspace(s1,s2,100);
ct12=linspace(ct1,ct2,100);
p12=linspace(p1,p2,100);

p13=0.5*(p1+p3);
p23=0.5*(p2+p3);

s3_=s3*ones(size(p12));
ct3_=ct3*ones(size(p12));
pref=0.5*(p12+p3);
F=gsw_rho(s12,ct12,pref)-gsw_rho(s3_,ct3_,pref);

sz=1.3*[10 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)])
plot(F,-p12)
hold on
plot([0 0],get(gca,'ylim')) 
hold off
ylabel('z [m]')
xlabel('F')
print('-dpdf','-r200','stab1_F.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_=linspace(33,38,50);
ct_=linspace(2,17,50);
[s_,ct_]=meshgrid(s_,ct_);

rho1=gsw_rho(s_,ct_,p1);
rho2=gsw_rho(s_,ct_,p2);
rho13=gsw_rho(s_,ct_,p13);
rho23=gsw_rho(s_,ct_,p23);


sz=1.3*[10 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)])
%contour(s_,ct_,rho1)
%hold on
contour(s_,ct_,rho1,gsw_rho([s1 s1],[ct1 ct1],[p1 p1]),'color','k');
hold on
contour(s_,ct_,rho2,gsw_rho([s2 s2],[ct2 ct2],[p2 p2]),'color','k');
contour(s_,ct_,rho13,gsw_rho([s1 s1],[ct1 ct1],[p13 p13]),'color','r','linestyle','--');
contour(s_,ct_,rho23,gsw_rho([s2 s2],[ct2 ct2],[p23 p23]),'color','r','linestyle','--');
plot(s1,ct1,'r*')
plot(s2,ct2,'b*')
plot(s3,ct3,'k*')
plot(s12,ct12)
hold off
ylabel('CT [deg C]')
xlabel('SA')
legend('ref local ','ref local',['ref mid-point ',num2str(pref(1))],['ref mid-point ',num2str(pref(end))],['cast (upper) ',num2str(p1)],['cast (lower) ',num2str(p2)],['bottle ',num2str(p3)],'cast lin. interp.','location','northwest')
print('-dpdf','-r200','stab1.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Method             SA       CT       P')
[SAns,CTns,pns] = depth_ntp_jackett(s3,ct3,p3,s12',ct12',p12');
disp(['Jackett:           ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])
[SAns,CTns,pns] = depth_ntp_guillaume(s3,ct3,p3,s12',ct12',p12');
disp(['Guillaume:         ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])
[SAns,CTns,pns] = depth_ntp_jackett_fzero(s3,ct3,p3,s12',ct12',p12');
disp(['Jackett   (fzero): ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])
[SAns,CTns,pns] = depth_ntp_guillaume_fzero(s3,ct3,p3,s12',ct12',p12');
disp(['Guillaume (fzero): ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])
[SAns,CTns,pns] = depth_ntp_iter(s3,ct3,p3,s12',ct12',p12');
disp(['depth_ntp_iter:    ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])

