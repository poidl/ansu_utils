
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='../exp276/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v259_1=drho_rms_hist; 
v259_2=drhoxy_rms_hist;

fname='../exp277/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v260_1=drho_rms_hist; 
v260_2=drhoxy_rms_hist;

fname='../exp278/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v270_1=drho_rms_hist; 
v270_2=drhoxy_rms_hist;

fname='../exp279/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v271_1=drho_rms_hist; 
v271_2=drhoxy_rms_hist;

nit=size(v259_1,1);

sz=1.4*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[0:nit-1];
p1=semilogy(xax,v259_2,'LineWidth',2)
hold on
semilogy(xax,v259_2,'+','LineWidth',3)
p2=semilogy(xax,v259_1,'--')
semilogy(xax,v259_1,'o')

p3=semilogy(xax,v260_2,'r','LineWidth',2)
semilogy(xax,v260_2,'r+','LineWidth',3)
p4=semilogy(xax,v260_1,'r--')
semilogy(xax,v260_1,'ro')

p5=semilogy(xax,v270_2,'g','LineWidth',2)
semilogy(xax,v270_2,'g+','LineWidth',3)
p6=semilogy(xax,v270_1,'g--')
semilogy(xax,v270_1,'go')

p7=semilogy(xax,v271_2,'k','LineWidth',2)
semilogy(xax,v271_2,'k+','LineWidth',3)
p8=semilogy(xax,v271_1,'k--')
semilogy(xax,v271_1,'ko')

power=[-11:-2];
set(gca, 'YTick', (10*ones(1,1:length(power))).^power )
grid on
legend([p1 p3  p5  p7 ],{'\sigma_0, b=1','\sigma_0, b variable',...
    '\sigma_p, b=1', '\sigma_p, b variable'});
title('helicity free')
xlabel('Iterations')
ylabel('$\rm kg/m^{-3}$','interpreter','latex','fontsize',13)
axis tight
print('-dpdf','-r200',['figures/compare_drho_rms_new_start'])