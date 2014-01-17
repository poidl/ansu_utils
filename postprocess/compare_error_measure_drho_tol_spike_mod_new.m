
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

runs={'exp315','exp316','exp317','exp318'};

drho=[];
drhoxy=[];

for i=1:length(runs)
	fname=['../',runs{i},'/data/iteration_history.mat'];
    load(fname,'drho_rms_hist','drhoxy_rms_hist' );
    
    drho=[drho;drho_rms_hist'];
    drhoxy=[drhoxy;drhoxy_rms_hist'];
    


end

mini=min(drhoxy');

nexp=size(drho,1);
nit=size(drho,2);


sz=1.4*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 

xax=[0:nit-1];

p1=nan*ones(1,nexp);
p2=nan*ones(1,nexp);
cols={'b','r','g','k'};

for i=1:nexp
    if 1
        p1(i)=semilogy(xax,drhoxy(i,:)-mini(1),'color',cols{i},'LineWidth',2)
        hold on
        semilogy(xax,drhoxy(i,:)-mini(1),'color',cols{i},'linestyle','+','LineWidth',2)
    else
        p1(i)=semilogy(xax,drhoxy(i,:),'color',cols{i},'LineWidth',2)
        hold on
        semilogy(xax,drhoxy(i,:),'color',cols{i},'linestyle','+','LineWidth',2)
    end
    
    p2(i)=semilogy(xax,drho(i,:),'color',cols{i},'linestyle','--','LineWidth',2)
    hold on
    semilogy(xax,drho(i,:),'color',cols{i},'linestyle','+','LineWidth',2)
end

power=[-11:-2];
set(gca, 'YTick', (10*ones(1,1:length(power))).^power )
grid on
legend(p1,{'tol=1e-15','tol=1e-9','tol=1e-8','tol=1e-7'})
    title('helicity spike')
    
xlabel('Iterations')
ylabel('$\rm kg/m^{-3}$','interpreter','latex','fontsize',13)
axis tight
print('-dpdf','-r200',['figures/compare_drho_rms_tol_spike_mod_new'])