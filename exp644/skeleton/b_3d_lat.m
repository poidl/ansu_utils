clear all;
close all;

addpath(genpath('../../../../gsw_matlab_v3_02'))
%addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))
omega_user_input; % for zonally_periodic

runs=450:548;

fname=['../../../ansu_utils_old/exp',num2str(runs(1)),'/data/iteration_history.mat'];
load(fname,'pns_hist');

fname=['../../../ansu_utils_old/exp',num2str(runs(1)),'/data/stationindex.mat'];
load(fname); % istation

load skeleton.mat 

b_surf=nan*ones(size(pns_3d));

if 0
    cnt=0;
    for ii=runs
        if mod(ii,5)==0
            disp(['reading exp',num2str(ii)])
        end
        cnt=cnt+1;
        fname=['../../../ansu_utils_old/exp',num2str(ii),'/data/input_data.mat'];
        load(fname,'s','ct','p');

        pns=squeeze(pns_3d(cnt,:,:));
        sns=squeeze(sns_3d(cnt,:,:));
        ctns=squeeze(ctns_3d(cnt,:,:));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [zi,yi,xi]=size(s);
        n2ns=n2(s,ct,p,pns);

        sns(isnan(n2ns))=nan;
        ctns(isnan(n2ns))=nan;
        pns(isnan(n2ns))=nan;

        n2nsx=0.5*(circshift(n2ns, [0 -1])+n2ns);
        if ~zonally_periodic;
            n2nsx(:,xi) = nan;
        end
        n2nsy=0.5*(circshift(n2ns, [-1 0])+n2ns);
        n2nsy(yi,:) = nan;

%         ro=gsw_rho(sns,ctns,pns);
%         rox=0.5*(circshift(ro, [0 -1])+ro);
%         if ~zonally_periodic;
%             rox(:,xi) = nan;
%         end
%         roy=0.5*(circshift(ro, [-1 0])+ro);
%         roy(yi,:) = nan;

        [kx_ns,ky_ns]=rho_p_xy(s,ct,p,pns);

        fx=9.81^2*(1./n2nsx).*kx_ns;
        fy=9.81^2*(1./n2nsy).*ky_ns;

        iyeq=(~isnan(pns) & ~isnan(circshift(pns,[-1 0])));
        bad=(iyeq & isnan(ky_ns));
        pns(bad)=nan;
        sns(bad)=nan;
        ctns(bad)=nan;
    
        ixeq=(~isnan(pns) & ~isnan(circshift(pns,[0 -1])));
        bad=(ixeq & isnan(kx_ns));
        pns(bad)=nan;
        sns(bad)=nan;
        ctns(bad)=nan;                

        regions=find_regions(pns);
        
        % find independent regions -> a least-squares problem is solved for
        % each of these regions
        
%         ixeq=(~isnan(n2ns) & ~isnan(circshift(n2ns,[0 -1])));
%         iyeq=(~isnan(n2ns) & ~isnan(circshift(n2ns,[-1 0])));
%         i1=(ixeq & ~isnan(fx));
%         i2=(iyeq & ~isnan(fy));
%         
%         i1(~i1 & circshift(i1,[0 1]))=true;
%         i2(~i2 & circshift(i2,[1 0]))=true;
%         
%         good=nan*ones(size(i1));
%         good(i1 & i2)=1;
%         regions=find_regions(good);

        lnb=solve_lsqr(regions,fx,fy);
        bns=exp(lnb);

        setnan=true(size(sns));
        for iregion=1:length(regions)
            region=regions{iregion};
            if ismember(istation,region)
                setnan(region)=false;
                bns(setnan)=nan;
                b_surf(cnt,:,:)=bns;

                found=true;
            end
        end

        if found==false
            error('hoitaus')
        end
        
%         if cnt==39
%             keyboard
%         end

    end
    save('b_surf.mat','b_surf')
else
    load('b_surf.mat')
    fname=['../../../ansu_utils_old/exp',num2str(runs(1)),'/data/input_data.mat'];
    load(fname,'p');
    [zi,yi,xi]=size(p);
    b=nan*ones(size(p));
    for ii=1:xi
        for jj=1:yi
            ps=squeeze(pns_3d(:,jj,ii));
            bs=squeeze(b_surf(:,jj,ii));
%             if any(~isnan(bs))
%                 keyboard
%             end
            pp=squeeze(p(:,jj,ii));
            for kk=1:zi
                above=find(ps<=pp(kk),1,'last');
                below=find(ps>pp(kk),1,'first');
                if ~isempty(above) && ~isempty(below) && below==(above+1)
                    dl=(pp(kk)-ps(below))/(ps(above)-ps(below));
                    b(kk,jj,ii)=bs(above)+dl*(bs(below)-bs(above));
                end
            end   
        end
    end
    
    %b=reshape(b,[zi,yi,xi]);
    save('../data/b.mat','b') 
    ncname='b.nc'
    delete(ncname)        
    nccreate(ncname,'b',...
             'Dimensions', {'k', zi, 'j', yi, 'i', xi})
    ncwrite(ncname,'b', b);
    
end




