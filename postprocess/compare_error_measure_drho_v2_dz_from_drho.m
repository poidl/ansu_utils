
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='../exp259/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v259_1=drho_rms_hist; 
v259_2=drhoxy_rms_hist;

fname='../exp269/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v260_1=drho_rms_hist; 
v260_2=drhoxy_rms_hist;

fname='../exp270/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v270_1=drho_rms_hist; 
v270_2=drhoxy_rms_hist;

fname='../exp271/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v271_1=drho_rms_hist; 
v271_2=drhoxy_rms_hist;

nit=size(v259_1,1);

sz=1.4*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[0:nit-1];
p1=semilogy(xax,v259_2)
hold on
semilogy(xax,v259_2,'+','LineWidth',2)
p2=semilogy(xax,v259_1,'--')
semilogy(xax,v259_1,'+','LineWidth',2)

p3=semilogy(xax,v260_2,'r')
semilogy(xax,v260_2,'r+','LineWidth',2)
p4=semilogy(xax,v260_1,'r--')
semilogy(xax,v260_1,'r+','LineWidth',2)

p5=semilogy(xax,v270_2,'g')
semilogy(xax,v270_2,'g+','LineWidth',2)
p6=semilogy(xax,v270_1,'g--')
semilogy(xax,v270_1,'g+','LineWidth',2)

p7=semilogy(xax,v271_2,'k')
semilogy(xax,v271_2,'k+','LineWidth',2)
p8=semilogy(xax,v271_1,'k--')
semilogy(xax,v271_1,'k+','LineWidth',2)

power=[-11:-2];
set(gca, 'YTick', (10*ones(1,1:length(power))).^power )
grid on
legend([p1 p2 p3 p4 p5 p6 p7 p8],{'no modification (using p_{ans}): rms gradient magnitude','no modification (using p_{ans}): rms potential','no b using p_{mid}: rms gradient magnitude',...
    'no b using p_{mid}: rms potential','b using p_{ans}: rms gradient magnitude','b using p_{ans}: rms potential', 'b using p_{mid}: rms gradient magnitude','b using p_{mid}: rms potential'});
title('helicity free')
xlabel('Iterations')
ylabel('$\rm kg/m^{-3}$','interpreter','latex','fontsize',13)
axis tight
print('-dpdf','-r200',['figures/compare_drho_rms_v2_dz_from_drho'])