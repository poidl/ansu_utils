% plotting script
clear all;
close all;

load('ansu_hist.mat',  'depth_change_e_i_hist');
load('/home/z3439823/mymatlab/omega/data_stefan/gk_ak_gamma_subs4x4.mat', 'lats','longs')

lat=squeeze(lats(1,:,1));
lon=squeeze(longs(1,1,:));

vv=depth_change_e_i_hist; % variable to plot
nit=size(vv,1);

nfig=3; % number of figures (pages)
ncols=2; % number of columns
nrows=3; % number of rows
nsp=ncols*nrows; % max. number of subplots per figure
%nsp=1
iit=1; % index of iteration

spwidth=0.4; % subplot width
spheight=0.2; % subplot height
wscol= 0.05; % white space between columns
wsrow=0.1; % white space between rows
leftmarg=(1-(ncols*spwidth+(ncols-1)*wscol))*0.5; % left and right margin
topmarg=(1-(nrows*spheight+(nrows-1)*wsrow))*0.5; % top and bottom margin

cmax=log10(max(abs(vv(:))));
%cmax=-5
cmin=-8
for ifig=1:nfig
    figure('PaperType','A4','PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
    isp=1; % index of subplot
    while (isp<=nsp) && ( ((ifig-1)*nsp+isp)<=nit)
        
        irow=ceil(isp/ncols); % row index
        icol=isp-(irow-1)*ncols; % column index
        disp(irow)
        disp(icol)
        left=leftmarg+(icol-1)*(spwidth+wscol); % current subplot position
        bottom=1-topmarg-irow*spheight-(irow-1)*wsrow; % current subplot position
      
        subplot('position',[left,bottom,spwidth,spheight])
        
        vp=squeeze(vv(iit,:,:));
        if 1
            tmp=vp
            tmp(vp<=0)=nan;
            pos=log10(   tmp );
            h=imagesc(lon,lat,pos);
            set(h,'alphadata',~isnan(pos)) % white nans
            set(gca,'YDir','normal')
            cm=colormap(flipud(colormap('hot'))) ;
            caxis([cmin cmax])
            cb=colorbar('location','southoutside','position',[left,bottom-0.4*wsrow,spwidth,0.1*wsrow],'XTickLabel',[])
            cbfreeze(cb)
            freezeColors
            hold on
        end
        if 1
            tmp=vp
            tmp(vp>=0)=nan;
            neg=log10(  -tmp );
            h=imagesc(lon,lat,neg);
            set(h,'alphadata',~isnan(neg)) % white nans
            set(gca,'YDir','normal')
            cm=colormap(fliplr(flipud(colormap('hot')))) ;
            caxis([cmin cmax])
            cb=colorbar('location','southoutside','position',[left,bottom-0.6*wsrow,spwidth,0.1*wsrow])
            xlabel(cb,'log10(depth change [m]); Red: pos. Blue: neg.')
        end

        

        iit=iit+1; isp=isp+1;
    end
    %colorbar('location','southoutside')
    print('-dpdf','-r500',['/home/z3439823/mymatlab/omega/figures/depth_change_e_i_hist/fig',num2str(ifig,'%02i')])
end
