
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='../exp324/data/iteration_history.mat';
load(fname,'drho_rms_hist','epsilon_rms_hist');
amp1=drho_rms_hist; 
eps1=epsilon_rms_hist;

e_final=min(eps1);



nit=size(amp1,1);

sz=1.4*[13 5];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[0:nit-1];

xax=xax(4:end);
eps1=eps1(4:end);
p1=semilogy(xax,eps1,'LineWidth',1)
hold on
semilogy(xax,eps1,'o','LineWidth',1)
grid on
xlabel('Iterations')
ylabel('$$\epsilon_{rms}^i\,\, \rm [kg/m^{-2}]$$','interpreter','latex','fontsize',13)
%text(-0.17,1.1,['a)'],'units','normalized','fontsize',12)

axis tight
print('-dpdf','-r200',['figures/paper_error_zoom'])
