
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

run1=571;
run2=574;

fname=['../exp',num2str(run1),'/data/iteration_history.mat'];
load(fname,'drho_rms_hist','epsilon_rms_hist','pns_hist');
amp1=drho_rms_hist; 
eps1=epsilon_rms_hist;
pr1=pns_hist;

fname=['../exp',num2str(run2),'/data/iteration_history.mat'];
load(fname,'drho_rms_hist','epsilon_rms_hist','pns_hist');
amp2=drho_rms_hist; 
eps2=epsilon_rms_hist;
pr2=pns_hist;

iend=10;
amp1=amp1(1:iend,:,:);
eps1=eps1(1:iend,:,:);
pr1=pr1(1:iend+1,:,:);
amp2=amp2(1:iend,:,:);
eps2=eps2(1:iend,:,:);
pr2=pr2(1:iend+1,:,:);

nit=size(amp1,1);

sz=1.2*[16 12];
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
xlim([1 4])
xlabel('Iterations')
ylabel('$$\epsilon_{rms}^i\,\, \rm [kg/m^{-2}]$$','interpreter','latex','fontsize',13)
ax2=axes('position',p);
set(ax2,'Visible','off')
text(-0.32,1.07,['a)'],'units','normalized','fontsize',12,'Parent',ax2)

pnew=p;
hfac=0.43;
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
meps2=min(eps2(:));
ylim([meps2,eps1(end)+2*10e-15])
ylim([meps2,eps1(end)+abs(eps1(end)-meps1)])

undershoot=abs(eps1(end)-meps1)
undershoot*1852*60*4

grid on
xlabel('Iterations')
ylabel('$$\epsilon_{rms}^i\,\, \rm [kg/m^{-2}]$$','interpreter','latex','fontsize',13)
text(-0.32,0.5,['c)'],'units','normalized','fontsize',12,'Parent',ax2)
xlim([1 nit])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sub=subplot(1,2,2);
p = get(sub, 'position');
hsub=gca;
p2=semilogy(xax,amp1,'');
hold on
semilogy(xax,amp1,'o')
p3=semilogy(xax,amp2,'r','LineWidth',1);
semilogy(xax,amp2,'ro','LineWidth',1)
xlim([1 nit])

grid on
xlabel('Iterations')
ylabel('$$\Phi_{rms}^i\,\, \rm [kg/m^{-3}]$$','interpreter','latex','fontsize',13)

ax2=axes('position',p);
set(ax2,'Visible','off')
%text(-0.32,1.07,['a)'],'units','normalized','fontsize',12,'Parent',ax2)
text(-0.22,1.07,['b)'],'units','normalized','fontsize',12,'Parent',ax2)

pnew=p;
hfac=0.43;
pnew(2) = p(2)+(1-hfac)*p(4);
pnew(4) = hfac*p(4);
set(sub, 'position', pnew);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax=axes('position',[p(1),p(2),p(3),hfac*p(4)]);

final=repmat(pr1(end,:,:),[size(pr1,1),1,1]);

pr1=abs(pr1-final);
pr1=max(pr1,[],3);
pr1=max(pr1,[],2);

pr2=abs(pr2-final);
pr2=max(pr2,[],3);
pr2=max(pr2,[],2);

pr1=pr1(1:end-1);
pr2=pr2(1:end-1);

p1=semilogy(xax,pr1,'LineWidth',1);
hold on
semilogy(xax,pr1,'o','LineWidth',1);
p4=semilogy(xax,pr2,'r');
semilogy(xax,pr2,'ro')

grid on
xlabel('Iterations')
%ylabel('$$\epsilon_{rms}^i\,\, \rm [kg/m^{-2}]$$','interpreter','latex','fontsize',13)
ylabel('$$\max_A{|p''^{\,i}|}\,\, \rm [db]$$','interpreter','latex','fontsize',13) % ,'fontsize',13)
text(-0.22,0.5,['d)'],'units','normalized','fontsize',12,'Parent',ax2)
xlim([1 iend])


print('-dpdf','-r200',['figures/paper_error_new_b',num2str(run1),'_',num2str(run2),'_mod'])