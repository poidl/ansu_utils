% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

dir={'041','043','042'};
dir2={'050','051','052'};

varname= 'eps_ss';

%figure('Visible','off')
nit=15
vv=nan(nit,length(dir));
xax=[0:size(vv)-1];
inds=4:14;
xaxp=xax(inds)
if 1;
    for ii=1:length(dir2);
        fname=['../exp',dir{ii},'/data/ansu_hist.mat'];
        load(fname, 'ithist');
        gr=ithist.(varname);
        size(gr)
        vv(:,ii)=ithist.(varname); % variable to plot
    end
    vvp=vv(inds,:)
    semilogy(xaxp,vvp,'o','MarkerFaceColor',[1. 1. 1.],'MarkerSize',10)
end
hold on
for ii=1:length(dir);
    fname=['../exp',dir{ii},'/data/ansu_hist.mat'];
    load(fname, 'ithist');
    gr=ithist.(varname);
    size(gr)
    vv(:,ii)=ithist.(varname); % variable to plot
end

vvp=vv(inds,:)
semilogy(xaxp,vvp)
hold on 
semilogy(xaxp,vvp,'+','LineWidth',2)

%legend(dir)


axis tight
    


xlabel('Iterations')
ylabel('$\sqrt{\overline{\epsilon^2}}\quad [m^{-1}]$','interpreter','latex','fontsize',13)
grid on

%print('-dpng','-r200',['../figures/',varname])
print('-dpng','-r200',['../figures/',varname,'_closeup'])
