% plotting script
clear all;
close all;

fname='data/iteration_history.mat';
varname= 'pns_hist';
varname1= 'ex_hist';
varname2= 'ey_hist';
load(fname, varname);
load(fname, varname1);
load(fname, varname2);
load('data/input_data.mat', 'lats','longs');

lat=squeeze(lats(1,:,1));
lon=squeeze(longs(1,1,:));

vv=pns_hist; % variables to plot
vv1=ex_hist; 
vv2=ey_hist;
nit=size(vv,1);

% averaging of ex/ey onto rho-points
%vv1=0.5*(vv1+circshift(vv1,[0,0,-1])) % assumes zonally periodic domain
%vv2=0.5*(vv2+circshift(vv2,[0,0,-1])) 

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

cbh=nan*ones(nsp,1); % colorbar handles
fac=nan*ones(nsp,1);
maxvecx=nan*ones(nsp,1);
maxvecy=nan*ones(nsp,1);

colorscale={'lin'};

for icolorscale=1:1
    iit=1; % index of iteration
    
    for ifig=1:nfig
        sz=2.0*[21 29.7];
        figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 
        set(gcf,'DefaultAxesFontSize', 15)
        set(gcf,'DefaultTextFontSize',15)

        isp=1; % index of subplot
        ip=1; % index of plot
        
        %cmp=colormap(hot(128));
        %cmp2=cmp(find(cmp(:,1)==1,1,'first')-15:end,:);
        
        while (isp<=nsp) && ( ip <=nit) && (iit<=size(vv,1))
            ip=((ifig-1)*nsp+isp);
            irow=ceil(isp/ncols); % row index
            icol=isp-(irow-1)*ncols; % column index
            left=leftmarg+(icol-1)*(spwidth+wscol); % current subplot position
            bottom=1-topmarg-irow*spheight-(irow-1)*wsrow; % current subplot position

            subplot('position',[left,bottom,spwidth,spheight])
            vp1=squeeze(vv1(iit,:,:));
            vp2=squeeze(vv2(iit,:,:));
            
            maxvecx(isp)=max(abs([vp1(:)]));
            maxvecy(isp)=max(abs([vp2(:)]));
            
            x=squeeze(longs(1,:,:));
            y=squeeze(lats(1,:,:));
            xx=0.5*(x+circshift(x,[0,-1]));
            yy=0.5*(y+circshift(y,[-1,0]));
            
            
            vp=squeeze(vv(iit,:,:));
            
            if colorscale{icolorscale}=='lin'
                tag='lin';
                
                h=imagesc(lon,lat,vp);
                hold on
                quiver(xx,y,vp1,0*vp1,'color','w','AutoScaleFactor',1.0,'LineWidth',3) % x-component
                quiver(x,yy,0*vp2,vp2,'color','k','AutoScaleFactor',1.0,'LineWidth',3) % y-component
                cbh(isp)=colorbar('Location','SouthOutside','position',[left,bottom-0.4*wsrow,spwidth,0.1*wsrow]);
                 if isp==6 % 
                     for ii=1:nsp
%                         set(cbh(ii),'XTick',-1:0.25:1)
%                         set(cbh(ii),'XTickLabel',num2str(fac(ii)*str2num(get(cbh(ii),'XTickLabel')),'%2.1e'));
%                         %set(cbh(ii),'XTickLabel',num2str(cmax2*str2num(get(cbh(ii),'XTickLabel')),'%2.1e'));
                         xlabel(cbh(ii), ['p_{ns} [db]   Max. length x-comp.:',num2str(maxvecx(ii)) ,', y-comp.: ',num2str(maxvecy(ii))], 'fontsize',18)
                     end
                 end

            end

            set(h,'alphadata',~isnan(vp)) % white nans
            set(gca,'YDir','normal')
            %load ('coast_data.mat');
            %plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);
            title(['Iteration ',num2str(ip-1)])
   


            iit=iit+1; isp=isp+1;
        end

        print('-dpng','-r200',['figures/epsilon_vectors_',num2str(ifig,'%02i')])
    end
end
