% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='data/iteration_history.mat';
varname= 'derr_rms_hist';
load(fname, varname);

vv1=derr_rms_hist; % variable to plot
nit=size(vv1,1);
vv1=nan*vv1;
varname= 'derrxy_rms_hist';
load(fname, varname);
vv2=derrxy_rms_hist; % variable to plot

figure('Visible','off')
xax=[0:size(vv1)-1];
semilogy(xax,vv1)
hold on
semilogy(xax,vv1,'+','LineWidth',2)
semilogy(xax,vv2,'--')
semilogy(xax,vv2,'+','LineWidth',2)
xlabel('Iterations')
ylabel('$$\sqrt{\overline{derr^2}}\quad \rm [kg/m^3]$$','interpreter','latex','fontsize',13)
grid on
print('-dpng','-r200',['figures/derr_derrxy_rms'])
