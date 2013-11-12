
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

fname='../exp259/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v259_1=drho_rms_hist; 
v259_2=drhoxy_rms_hist;

fname='../exp260/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v260_1=drho_rms_hist; 
v260_2=drhoxy_rms_hist;

fname='../exp261/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v261_1=drho_rms_hist; 
v261_2=drhoxy_rms_hist;

fname='../exp262/data/iteration_history.mat';
load(fname,'drho_rms_hist','drhoxy_rms_hist' );
v262_1=drho_rms_hist; 
v262_2=drhoxy_rms_hist;

nit=size(v259_1,1);

sz=1.4*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[0:nit-1];
% p1=semilogy(xax,v259_2)
% hold on
% semilogy(xax,v259_2,'+','LineWidth',2)
% p2=semilogy(xax,v259_1,'--')
% semilogy(xax,v259_1,'+','LineWidth',2)
% 
% p3=semilogy(xax,v260_2,'r')
% semilogy(xax,v260_2,'r+','LineWidth',2)
% p4=semilogy(xax,v260_1,'r--')
% semilogy(xax,v260_1,'r+','LineWidth',2)

p1=semilogy(xax,v261_2)
hold on
semilogy(xax,v261_2,'+','LineWidth',2)
p2=semilogy(xax,v261_1,'--')
semilogy(xax,v261_1,'+','LineWidth',2)
p5=semilogy(xax,v259_2,'o','MarkerSize',8)

p3=semilogy(xax,v262_2,'r')
semilogy(xax,v262_2,'r+','LineWidth',2)
p4=semilogy(xax,v262_1,'r--')
semilogy(xax,v262_1,'r+','LineWidth',2)
p6=semilogy(xax,v260_2,'ro','MarkerSize',8)

power=[-11:-2];
set(gca, 'YTick', (10*ones(1,1:length(power))).^power )
grid on
%legend([p1 p2 p3 p4],{'no modification: rms gradient magnitude','no modification: rms potential','using b: rms gradient magnitude','using b: rms potential'});
legend([p1 p2 p3 p4 p5 p6],{'no modification: rms gradient magnitude','no modification: rms potential', ...
                                'using b: rms gradient magnitude','using b: rms potential', 'helicity free rms potential','helicity free rms potential'},...
                                'location',[0.6 0.5 0.2 0.2]);
title('helicity spike')
xlabel('Iterations')
ylabel('$\rm kg/m^{-3}$','interpreter','latex','fontsize',13)
axis tight
print('-dpdf','-r200',['figures/compare_drho_rms_v2b'])