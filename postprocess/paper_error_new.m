
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='../exp568/data/iteration_history.mat';
load(fname,'drho_rms_hist','epsilon_rms_hist');
amp1=drho_rms_hist; 
eps1=epsilon_rms_hist;

fname='../exp569/data/iteration_history.mat';
load(fname,'drho_rms_hist','epsilon_rms_hist');
amp2=drho_rms_hist; 
eps2=epsilon_rms_hist;


nit=size(amp1,1);

sz=1.4*[16 8];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[1:nit];

sub=subplot(1,2,1);
p = get(sub, 'position');
hsub=gca;
p1=semilogy(xax,eps1,'LineWidth',1);
hold on
semilogy(xax,eps1,'o','LineWidth',1)
p4=semilogy(xax,eps2,'r');
semilogy(xax,eps2,'ro')
grid on
xlim([1 nit])
xlabel('Iterations')
ylabel('$$\epsilon_{rms}^i\,\, \rm [kg/m^{-2}]$$','interpreter','latex','fontsize',13)
ax2=axes('position',p);
set(ax2,'Visible','off')
text(-0.32,1.07,['a)'],'units','normalized','fontsize',12,'Parent',ax2)

pnew=p;
hfac=0.4;
pnew(2) = p(2)+(1-hfac)*p(4);
pnew(4) = hfac*p(4);
set(sub, 'position', pnew);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax=axes('position',[p(1),p(2),p(3),hfac*p(4)]);
p1=semilogy(xax,eps1,'LineWidth',1);
hold on
semilogy(xax,eps1,'o','LineWidth',1);
p4=semilogy(xax,eps2,'r');
semilogy(xax,eps2,'ro')

meps1=min(eps1(:));
ylim([meps1,eps1(end)+2*10e-15])
ylim([meps1,eps1(end)+abs(eps1(end)-meps1)])

undershoot=abs(eps1(end)-meps1)
undershoot*1852*60*4

grid on
xlabel('Iterations')
ylabel('$$\epsilon_{rms}^i\,\, \rm [kg/m^{-2}]$$','interpreter','latex','fontsize',13)
text(-0.32,0.5,['c)'],'units','normalized','fontsize',12,'Parent',ax2)
xlim([1 nit])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)
p2=semilogy(xax,amp1,'');
hold on
semilogy(xax,amp1,'o')
p3=semilogy(xax,amp2,'r','LineWidth',1);
semilogy(xax,amp2,'ro','LineWidth',1)
text(-0.22,1.07,['b)'],'units','normalized','fontsize',12)
xlim([1 nit])
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
print('-dpdf','-r200',['figures/paper_error_new'])
