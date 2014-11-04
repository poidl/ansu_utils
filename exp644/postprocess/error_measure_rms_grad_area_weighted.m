% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='data/iteration_history.mat';
varname= 'err_grad_rms_area_weighted_hist';
load(fname, varname);

vv1=err_grad_rms_area_weighted_hist; % variable to plot

figure('Visible','off')
xax=[0:size(vv1)-1];
semilogy(xax,vv1)
hold on
semilogy(xax,vv1,'+','LineWidth',2)
xlabel('Iterations')
%ylabel('$$\sqrt{\overline{derr^2}}\quad \rm [kg/m^3]$$','interpreter','latex','fontsize',13)
ylabel('area weighted rms gradient')
grid on
print('-dpng','-r200',['figures/rms_gradient_area_weighted'])