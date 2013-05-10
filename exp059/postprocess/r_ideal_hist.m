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

it=2
v=r_idea(it,:,:)
hist(v(:),50)
ylim([0 5])
%[nelements,xcenters]=hist(v(:),50);
%bar(xcenters,log10(nelements))
xlabel('$\tilde{r}$','interpreter','latex','fontsize',15)
ylabel('frequency','fontsize',15)


print('-dpng','-r200',['figures/r_ideal_hist'])

