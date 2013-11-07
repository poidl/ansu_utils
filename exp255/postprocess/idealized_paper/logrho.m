clear all
close all

load iteration_history

sz=1.0*[10 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 
%set(gcf,'DefaultAxesFontSize', 15)
%set(gcf,'DefaultTextFontSize',15)

drho=drho_hist(:,2,1);


ipos=drho>0;
ineg=drho<0;

%drho_p=nan*drho;
%drho_p(ipos)=drho(ipos)

%drho_n=nan*drho;
%drho_n(ineg)=drho(ineg)

logr=nan*drho;
logr(ipos)=log10(drho(ipos));
logr(ineg)=-log10(-drho(ineg));
%logr(logr==-Inf)=nan;

plot(logr,'*')
xlabel('Iteration Nr')
ylabel('$\pm log10( \pm \tilde{\Delta}\rho ) \rm\,[kg/m^3]$','interpreter','latex')
xlim([0.8 8.5])
set(gca,'XTick',[1:9])
grid on
hold on


%%%%%%%%%% first term in Eq. 18
sns=squeeze( sns_hist(:,2,1) );
ctns=squeeze( ctns_hist(:,2,1) );
pns=squeeze( ctns_hist(:,2,1) );

Delta_plus1=nan*ones(1,length(sns)-1);

for it=1:9;
    pmid1=0.5*(pns_hist(it,1,1)+pns_hist(it,2,1));
    %pmid2=0.5*(pns_hist(it+1,1,1)+pns_hist(it+1,2,1));

    rBpi=gsw_rho(sns(it),ctns(it),pmid1);
    rCpi=gsw_rho(sns(it+1),ctns(it+1),pmid1);

    rBpi_=gsw_rho(sns(it),ctns(it),pmid1+1);
    rCpi_=gsw_rho(sns(it+1),ctns(it+1),pmid1+1);

    Delta_plus1(it)=((rBpi_-rBpi)-(rCpi_-rCpi))*500;
end

%Delta_plus1=[nan, Delta_plus1]
Delta_plus1=Delta_plus1(1:length(logr))
%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(log10(Delta_plus1),'ro')



print('-dpdf','-r200','figures/logrho.pdf')