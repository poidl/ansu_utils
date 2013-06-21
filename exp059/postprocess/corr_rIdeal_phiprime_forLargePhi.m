clear all
close all

fname='data/iteration_history.mat';
varname= 'phiprime_e_hist';
load(fname, varname);
load('data/input_data.mat', 'lats','longs'); 

pp=phiprime_e_hist;
r=1;

summ=0;

pp=-pp;

iit=1;
cp_iit=prod(1+r*pp(iit:end,:,:), 1);
pp_iit=pp(iit,:,:);
r_idea_iit=(cp_iit-1)./ pp_iit;

pp_iit=pp(iit,:,:);

sig=std(pp_iit(~isnan(pp_iit)));

keep=(abs(pp_iit)>=sig);
pp_iitlarge=pp_iit(keep);
r_idea_iit_large=r_idea_iit(keep);

lats=lats(1,:,:);
lat=lats(keep);


nans=isnan(pp_iitlarge(:))|isnan(r_idea_iit_large(:));
c_all=polyfit(pp_iitlarge(~nans),r_idea_iit_large(~nans),1); % returns [slope, intercept]
scatter(pp_iitlarge,r_idea_iit_large,22,lat,'fill')
axis tight
hold on
x=linspace(min(pp_iitlarge()), max(pp_iitlarge()));
plot(x,c_all(1)*x+c_all(2))
text(0.1,1.0,['Slope: ',num2str(c_all(1),'%10.1e'),' Intercept: ', num2str(c_all(2))],'fontsize',15','units','normalized')

plot(get(gca,'xlim'), [1 1],'k','LineWidth',1)
plot( [0 0],get(gca,'ylim'),'k','LineWidth',1)


caxis([-90 90])
cb=colorbar()
ylabel(['$\tilde{r_{1}}$'],'interpreter','latex','fontsize',20)
xlabel('$\Phi^{''1}$','interpreter','latex','fontsize',15)
ylabel(cb,'Latitude (deg)')



print('-dpng','-r200',['figures/corr_rIdeal_phiprime_forLargePhi'])
