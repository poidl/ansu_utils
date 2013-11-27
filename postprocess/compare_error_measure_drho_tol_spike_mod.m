
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='../exp301/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v259_1=drho_rms_hist; 
v259_2=drhoxy_rms_hist;

e_final=min(v259_2);

fname='../exp302/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v260_1=drho_rms_hist; 
v260_2=drhoxy_rms_hist;

fname='../exp303/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v270_1=drho_rms_hist; 
v270_2=drhoxy_rms_hist;

fname='../exp304/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v271_1=drho_rms_hist; 
v271_2=drhoxy_rms_hist;

nit=size(v259_1,1);

sz=1.4*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[0:nit-1];
p1=semilogy(xax,v259_2-e_final,'LineWidth',2)
hold on
semilogy(xax,v259_2-e_final,'+','LineWidth',2)
p2=semilogy(xax,v259_1,'--')
semilogy(xax,v259_1,'o')

% p3=semilogy(xax,v260_2-e_final,'r','LineWidth',2)
% semilogy(xax,v260_2-e_final,'r+','LineWidth',2)
% p4=semilogy(xax,v260_1,'r--')
% semilogy(xax,v260_1,'ro')
% 
% p5=semilogy(xax,v270_2-e_final,'g','LineWidth',2)
% semilogy(xax,v270_2-e_final,'g+','LineWidth',2)
% p6=semilogy(xax,v270_1,'g--')
% semilogy(xax,v270_1,'go')
% 
% p7=semilogy(xax,v271_2-e_final,'k','LineWidth',2)
% semilogy(xax,v271_2-e_final,'k+','LineWidth',2)
% p8=semilogy(xax,v271_1,'k--')
% semilogy(xax,v271_1,'ko')

power=[-11:-2];
set(gca, 'YTick', (10*ones(1,1:length(power))).^power )
grid on
    title('helicity spike')
xlabel('Iterations')
ylabel('$\rm kg/m^{-3}$','interpreter','latex','fontsize',13)
axis tight
print('-dpdf','-r200',['figures/compare_drho_rms_tol_spike_mod'])
