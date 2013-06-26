clear all
close all

fname='data/iteration_history.mat';
varname= 'phiprime_e_hist';
load(fname, varname);
load('data/input_data.mat', 'lats','longs'); 

pp=phiprime_e_hist;
r=1;

sum=0;

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

pos=pp1_>0;
neg=pp1_<0;
semilogx(pp1_(pos),r1_21(pos),'b.','MarkerSize',12);
hold on
semilogx(-pp1_(neg),r1_21(neg),'r.','MarkerSize',12);

plot(get(gca,'xlim'), [1 1],'k--','LineWidth',1)


ll=logspace(-7.8,-3,100);
plot(ll,1+1e-6./(ll.^0.99));
plot(ll,1+1e-6./(-ll.^0.99),'r');


ylabel('$\tilde{r}$','interpreter','latex','fontsize',15)
xlabel('$\Phi^{''1}$','interpreter','latex','fontsize',15)


print('-dpng','-r200',['figures/corr_rIdeal_phiprime'])
