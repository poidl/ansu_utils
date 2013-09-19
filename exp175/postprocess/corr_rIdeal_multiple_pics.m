% plotting script
clear all;
close all;

fname='data/iteration_history.mat';
varname= 'phiprime_e_hist';
load(fname, varname);
load('data/input_data.mat', 'lats','longs'); 

pp=phiprime_e_hist;
r=1;

summ=0;

pp=-pp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nit=size(pp,1);

nfig=1; % number of figures (pages)
ncols=2; % number of columns
nrows=3; % number of rows
nsp=ncols*nrows; % max. number of subplots per figure
%nsp=1

spwidth=0.4; % subplot width
spheight=0.2; % subplot height
wscol= 0.08; % white space between columns
wsrow=0.1; % white space between rows
leftmarg=(1-(ncols*spwidth+(ncols-1)*wscol))*0.5; % left and right margin
topmarg=(1-(nrows*spheight+(nrows-1)*wsrow))*0.5; % top and bottom margin


sp=nan*ones(nsp,1); % colorbar handles
sig=nan*ones(nsp,1); 

colorscale={'lin'};
txt={'a)','b)','c)','d)','e)','f)'};
for icolorscale=1:1
    iit=1; % index of iteration
    
    for ifig=1:nfig
        sz=2.0*[21 29.7];
        figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off')
        %figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)])
        set(gcf,'DefaultAxesFontSize', 15)
        set(gcf,'DefaultTextFontSize',15)

        isp=1; % index of subplot
        ip=1; % index of plot
        
        while (isp<=nsp) && ( ip <=nit) && (iit<=nit)
            ip=((ifig-1)*nsp+isp);
            irow=ceil(isp/ncols); % row index
            icol=isp-(irow-1)*ncols; % column index
            left=leftmarg+(icol-1)*(spwidth+wscol); % current subplot position
            bottom=1-topmarg-irow*spheight-(irow-1)*wsrow; % current subplot position

            sp(isp)=subplot('position',[left,bottom,spwidth,spheight]);
            
            cp_iit=prod(1+r*pp(iit:end,:,:), 1);
            pp_iit=pp(iit,:,:);
            r_idea_iit=(cp_iit-1)./ pp_iit;

            pp_iit=pp(iit,:,:);

            sig(isp)=std(pp_iit(~isnan(pp_iit)));

            keep=(abs(pp_iit)>=sig(isp));
            pp_iitlarge=pp_iit(keep);
            r_idea_iit_large=r_idea_iit(keep);

            lats=lats(1,:,:);
            lat=lats(keep);


            nans=isnan(pp_iitlarge(:))|isnan(r_idea_iit_large(:));

            c_all=polyfit(pp_iitlarge(~nans),r_idea_iit_large(~nans),1); % returns [slope, intercept]
            %scatter(pp_iitlarge,r_idea_iit_large,22,lat,'fill')
            
            scatter(pp_iit(:),r_idea_iit(:),22,reshape(lats(1,:,:),[],1),'fill')
            axis tight
            hold on
            x=linspace(min(pp_iitlarge()), max(pp_iitlarge()));
            plot(x,c_all(1)*x+c_all(2))
            text(0.1,1.0,['Slope: ',num2str(c_all(1),'%10.1e'),' Intercept: ', num2str(c_all(2))],'fontsize',15','units','normalized')
            if isp==6 % 
                for ii=1:nsp
                    plot(sp(ii),get(sp(ii),'xlim'), [1 1],'k','LineWidth',1)
                    plot(sp(ii), [0 0],get(sp(ii),'ylim'),'k','LineWidth',1)
                    plot(sp(ii), [-sig(ii), -sig(ii)],get(sp(ii),'ylim'),'k--','LineWidth',1)
                    plot(sp(ii), [sig(ii), sig(ii)],get(sp(ii),'ylim'),'k--','LineWidth',1)
                end
            end

    
            cb=colorbar();
            caxis(90*[-1 1])
            ylabel(['$\tilde{r_{',num2str(iit),'}}$'],'interpreter','latex','fontsize',20)
            xlabel(['$\Phi^{''',num2str(iit),'}$'],'interpreter','latex','fontsize',20)
            ylabel(cb,'Latitude (deg)')
            title(['Iteration ',num2str(ip)])
            text(-20,95,txt(isp),'fontsize',18')

            iit=iit+1; isp=isp+1;
        end
        
        print('-dpdf','-r200',['figures/corr_rIdeal_multiple_pics_',num2str(ifig,'%02i')])
    end
end



