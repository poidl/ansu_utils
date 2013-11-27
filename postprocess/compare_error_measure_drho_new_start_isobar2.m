
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='../exp307/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
%v259_1=drho_rms_hist; 
l1=drhoxy_rms_hist;

fname='../exp308/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
%v260_1=drho_rms_hist; 
l2=drhoxy_rms_hist;

fname='../exp309/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
%v270_1=drho_rms_hist; 
l3=drhoxy_rms_hist;

fname='../exp310/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
%v271_1=drho_rms_hist; 
l4=drhoxy_rms_hist;


fname='../exp311/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
%v311_1=drho_rms_hist; 
l5=drhoxy_rms_hist;

fname='../exp312/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
%v312_1=drho_rms_hist; 
l6=drhoxy_rms_hist;

fname='../exp313/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
%v313_1=drho_rms_hist; 
l7=drhoxy_rms_hist;

fname='../exp314/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
%v314_1=drho_rms_hist; 
l8=drhoxy_rms_hist;

nit=size(l1,1);

sz=1.4*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[0:nit-1];
p1=semilogy(xax,l1,'LineWidth',2)
hold on
semilogy(xax,l1,'+','LineWidth',3)


p3=semilogy(xax,l2,'--','LineWidth',2)
semilogy(xax,l2,'+','LineWidth',3)


p5=semilogy(xax,l3,'r','LineWidth',2)
semilogy(xax,l3,'r+','LineWidth',3)


p7=semilogy(xax,l4,'r--','LineWidth',2)
semilogy(xax,l4,'r+','LineWidth',3)


p8=semilogy(xax,l5,'g','LineWidth',2)
semilogy(xax,l5,'g+','LineWidth',3)


p9=semilogy(xax,l6,'g--','LineWidth',2)
semilogy(xax,l6,'g+','LineWidth',3)

p10=semilogy(xax,l7,'k','LineWidth',2)
semilogy(xax,l7,'k+','LineWidth',3)


p11=semilogy(xax,l8,'k--','LineWidth',2)
semilogy(xax,l8,'k+','LineWidth',3)

power=[-11:-2];
set(gca, 'YTick', (10*ones(1,1:length(power))).^power )
grid on
legend([p1 p3  p5  p7 p8 p9 p10 p11],{'p_{start}=1500, b=1','p_{start}=1500, b variable',...
    'p_{start}=1200, b=1', 'p_{start}=1200, b variable',...
    'p_{start}=900, b=1','p_{start}=900, b variable',...
    'p_{start}=800, b=1','p_{start}=800, b variable'});
title('helicity free')
xlabel('Iterations')
ylabel('$\rm kg/m^{-3}$','interpreter','latex','fontsize',13)
axis tight
print('-dpdf','-r200',['figures/compare_drho_rms_new_start_isobar2'])