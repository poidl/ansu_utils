
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='../exp568/data/iteration_history.mat';
load(fname,'pns_hist');
amp1=pns_hist;

fname='../exp569/data/iteration_history.mat';
load(fname,'pns_hist');
amp2=pns_hist; 

final=repmat(amp1(end,:,:),[size(amp1,1),1,1]);

amp1=abs(amp1-final);
amp1=max(amp1,[],3);
amp1=max(amp1,[],2);

amp2=abs(amp2-final);
amp2=max(amp2,[],3);
amp2=max(amp2,[],2);

nit=size(amp1,1);

sz=1.4*[16 8];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[1:nit];


p1=semilogy(xax,amp1,'LineWidth',1);
hold on
semilogy(xax,amp1,'o','LineWidth',1)
p4=semilogy(xax,amp2,'r');
semilogy(xax,amp2,'ro')
grid on
xlim([1 nit])
xlabel('Iterations')
%ylabel('$$\epsilon_{rms}^i\,\, \rm [kg/m^{-2}]$$','interpreter','latex','fontsize',13)


print('-dpdf','-r200',['figures/paper_error_pressure'])
