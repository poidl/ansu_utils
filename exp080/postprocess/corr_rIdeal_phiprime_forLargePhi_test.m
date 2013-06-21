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
cp=prod(1+r*pp, 1);

pp1=pp(1,:,:);
r1_21=(cp-1)./ pp1;
%cut=0.3;
%cutoff=r_idea<1-cut | r_idea>1+cut;
%r_idea(cutoff)=nan;

pp1_=pp(1,:,:);

%plot(pp1_(:),r1_21(:),'.')

sig=std(pp1_(~isnan(pp1_)));

keep=abs(pp1_)>=sig
pp1_large=pp1_(keep);
r1_21_large=r1_21(keep);

lats=lats(1,:,:);
lat=lats(keep);

longs=longs(1,:,:);
lon=longs(keep);

nans=isnan(pp1_large(:))|isnan(r1_21_large(:));
R = corrcoef(pp1_large(~nans),r1_21_large(~nans));
c_all=polyfit(pp1_large(~nans),r1_21_large(~nans),1) % returns [slope, intercept]
%plot(pp1_large,r1_21_large,'.')
scatter(pp1_large,r1_21_large,10,lat,'fill')
hold on
x=linspace(min(pp1_large()), max(pp1_large()));
plot(x,c_all(1)*x+c_all(2))
plot(x,1+750*x)
plot(get(gca,'xlim'), [1 1],'k','LineWidth',1)
plot( [0 0],get(gca,'ylim'),'k','LineWidth',1)
if 0
    cmp=colormap(hot(128));
    cmp2=cmp(find(cmp(:,1)==1,1,'first')-15:end,:);
    colormap([fliplr(cmp2);flipud(cmp2)]) ;
end
caxis([-90 90])
%caxis([0 360])
cb=colorbar()
ylabel('$\tilde{r}$','interpreter','latex','fontsize',15)
xlabel('$\Phi^{''1}$','interpreter','latex','fontsize',15)
ylabel(cb,'Latitude (deg)')


npos=sum(pp1_large(:)>0)
nneg=sum(pp1_large(:)<0)


print('-dpng','-r200',['figures/corr_rIdeal_phiprime_forLargePhi_test'])