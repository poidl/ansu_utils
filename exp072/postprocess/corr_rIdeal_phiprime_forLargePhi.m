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
cp=cumprod(1+r*pp, 1);

pp1=repmat(pp(1,:,:),[size(pp,1),1,1]);
r_idea=(cp-1)./ pp1;
%cut=0.3;
%cutoff=r_idea<1-cut | r_idea>1+cut;
%r_idea(cutoff)=nan;

%%%%%%%%%%%%%%%%%%%%
cp=cumprod(exp(r*pp), 1);
r_idea2=log(cp)./ pp1;

%%%%%%%%%%%%%%%%%%%%

r1_21=r_idea(21,:,:);
r2_21=r_idea2(21,:,:);


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
c_all=polyfit(pp1_large(~nans),r1_21_large(~nans),1); % p(1) slope, p(2) intercept
%plot(pp1_large,r1_21_large,'.')
scatter(pp1_large,r1_21_large,10,lat,'fill')
hold on
x=linspace(min(pp1_large()), max(pp1_large()));
%plot(x,c_all(1)*x+1)
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


% pos=pp1_>0;
% neg=pp1_<0;
% semilogx(pp1_(pos),r1_21(pos),'b.','MarkerSize',12);
% hold on
% semilogx(-pp1_(neg),r1_21(neg),'r.','MarkerSize',12);
% 
% plot(get(gca,'xlim'), [1 1],'k--','LineWidth',1)
% 
% 
% ll=logspace(-7.8,-3,100);
% plot(ll,1+1e-6./(ll.^0.99));
% plot(ll,1+1e-6./(-ll.^0.99),'r');
% 
% 
% ylabel('$\tilde{r}$','interpreter','latex','fontsize',15)
% xlabel('$\Phi^{''1}$','interpreter','latex','fontsize',15)


print('-dpng','-r200',['figures/corr_rIdeal_phiprime_forLargePhi'])
