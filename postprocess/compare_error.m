
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='../exp556/data/iteration_history.mat';
load(fname,'drho_rms_hist','epsilon_rms_hist');
amp1=drho_rms_hist; 
eps1=epsilon_rms_hist;

fname='../exp558/data/iteration_history.mat';
load(fname,'drho_rms_hist','epsilon_rms_hist');
amp2=drho_rms_hist; 
eps2=epsilon_rms_hist;


nit=size(amp1,1);

sz=1.4*[16 8];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[0:nit-1];

subplot(1,2,1)
p1=semilogy(xax,eps1,'LineWidth',1)
hold on
semilogy(xax,eps1,'o','LineWidth',1)
p4=semilogy(xax,eps2,'r')
semilogy(xax,eps2,'ro')
grid on
xlabel('Iterations')
ylabel('$$\epsilon_{rms}^i\,\, \rm [kg/m^{-2}]$$','interpreter','latex','fontsize',13)
text(-0.22,1.07,['a)'],'units','normalized','fontsize',12)



subplot(1,2,2)
p2=semilogy(xax,amp1,'')
hold on
semilogy(xax,amp1,'o')
p3=semilogy(xax,amp2,'r','LineWidth',1)
semilogy(xax,amp2,'ro','LineWidth',1)
text(-0.22,1.07,['b)'],'units','normalized','fontsize',12)
%set(gca,'ytick',10.^[-15:1:0])
%set(gca,'yticklabel'
%set(gca,'yminortick','off')


grid on
%set(gca,'yminorgrid','off')
%legend([p1  p3 p5 p7  ],{'tol=1e-15','tol=1e-11','tol=1e-8','tol=1e-5'})
%    title('full data set')
xlabel('Iterations')
ylabel('$$\Phi_{rms}^i\,\, \rm [kg/m^{-3}]$$','interpreter','latex','fontsize',13)
%axis tight
print('-dpdf','-r200',['figures/compare_error556_558'])
