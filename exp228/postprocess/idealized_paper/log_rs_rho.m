clear all
close all

load iteration_history

sz=1.0*[10 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 
%set(gcf,'DefaultAxesFontSize', 15)
%set(gcf,'DefaultTextFontSize',15)

drho=drho_hist(:,2,1);

rs_drho=sqrt(drho.^2);

log_rs=log10(rs_drho);

plot(log_rs,'*')
xlabel('Iteration Nr')
ylabel('$log10\left( \sqrt{\tilde{\Delta}\rho^2} \right) \rm\,[kg/m^3]$','interpreter','latex')
xlim([0.8 8.5])
ylim([-10.5 -2])
set(gca,'XTick',[1:9])
grid on
hold on


print('-dpdf','-r200','figures/log_rs_rho.pdf')