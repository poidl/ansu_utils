% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='data/iteration_history.mat';
varname= 'eps_rms_hist';
load(fname, varname);

vv=eps_rms_hist; % variable to plot
nit=size(vv,1);

figure('Visible','off')
xax=[0:size(vv)-1];
%xax=xax(3:end);
%vv=vv(3:end);
semilogy(xax,vv)
hold on
semilogy(xax,vv,'+','LineWidth',2)
xlabel('Iterations')
ylabel('$\sqrt{\overline{\epsilon^2}}\quad [m^{-1}]$','interpreter','latex','fontsize',13)
grid on
print('-dpng','-r200',['figures/eps_rms'])
