
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='../exp324/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
amp1=drho_rms_hist; 
eps1=drhoxy_rms_hist;

e_final=min(eps1);

fname='../exp325/data/iteration_history.mat';
load(fname,'drho_rms_hist','epsilon_rms_hist', 'drhoxy_rms_hist');
amp2=drho_rms_hist; 
eps2=epsilon_rms_hist;
drhoxy2=drhoxy_rms_hist;



nit=size(amp1,1);

sz=1.4*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[0:nit-1];

subplot(2,1,1)
p1=semilogy(xax,eps1,'LineWidth',1)
hold on
semilogy(xax,eps1,'o','LineWidth',1)
p4=semilogy(xax,amp2,'r--')
semilogy(xax,amp2,'ro')
grid on
xlabel('Iterations')
ylabel('$$\epsilon_{rms}^i \rm kg/m^{-3}$$','interpreter','latex','fontsize',13)

subplot(2,1,2)
p2=semilogy(xax,amp1,'r')
hold on
semilogy(xax,amp1,'ro')
p3=semilogy(xax,eps2-e_final,'--','LineWidth',1)
semilogy(xax,eps2-e_final,'o','LineWidth',1)



grid on
%legend([p1  p3 p5 p7  ],{'tol=1e-15','tol=1e-11','tol=1e-8','tol=1e-5'})
%    title('full data set')
xlabel('Iterations')
ylabel('$\rm kg/m^{-3}$','interpreter','latex','fontsize',13)
%axis tight
print('-dpdf','-r200',['figures/paper_error'])
