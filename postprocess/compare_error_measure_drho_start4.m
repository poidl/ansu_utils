% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='../exp251/data/iteration_history.mat';
varname= 'drho_rms_hist';
load(fname, varname);

vv1=drho_rms_hist; % variable to plot

fname='../exp256/data/iteration_history.mat';
varname= 'drho_rms_hist';
load(fname, varname);

vv2=drho_rms_hist; % variable to plot

nit=size(vv1,1);

sz=1.4*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[0:nit-1];
semilogy(xax,vv1)
hold on
semilogy(xax,vv1,'+','LineWidth',2)
semilogy(xax,vv2,'r')
semilogy(xax,vv2,'r+','LineWidth',2)
semilogy(xax(1:end-3),vv1(4:end),'b')
semilogy(xax(1:end-3),vv1(4:end),'b+','LineWidth',2)
power=[-11:-2];
set(gca, 'YTick', (10*ones(1,1:length(power))).^power )
grid on
xlabel('Iterations')
ylabel('$\sqrt{\overline{\Delta\rho^2}}\quad \rm [kg/m^{-3}]$','interpreter','latex','fontsize',13)
axis tight
print('-dpdf','-r200',['figures/compare_drho_rms_start4'])