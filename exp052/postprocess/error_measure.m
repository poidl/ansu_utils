% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='data/ansu_hist.mat';
varname= 'eps_ss';
load(fname, 'ithist');

vv=ithist.(varname); % variable to plot
nit=size(vv,1);

figure('Visible','off')
xax=[0:size(vv)-1];
semilogy(xax,vv)
hold on
semilogy(xax,vv,'+','LineWidth',2)
xlabel('Iterations')
ylabel('$\sqrt{\overline{\epsilon^2}}\quad [m^{-1}]$','interpreter','latex','fontsize',13)
grid on
print('-dpng','-r200',['figures/',varname])
