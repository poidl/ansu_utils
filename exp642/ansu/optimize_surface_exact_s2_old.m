function [sns_i,ctns_i,pns_i,rmsdrho,mdf,df_med] = optimize_surface_exact_s2(s,ct,p,sns,ctns,pns)

%           Optimize density surfaces to minimise the fictitious diapycnal diffusivity
%
%           Optimize density surface using an iterative method to minimise
%           fictitious diapycnal diffusivity according 'A new method of
%           forming approximately neutral surfaces, Klocker et al., Ocean
%           Science, 5, 155-172, 2009.'
%
%
% Input:    s           salinity
%           ct          conservative temperature
%           p           pressure
%           sns         salinity on initial density surface
%           ctns        conservative temperature on initial density
%                       surface
%           pns         pressure on initial density surface
%           e1t         grid length in meters in x-direction
%           e2t         grid length in meters in y-direction
%
% Output:   sns_i       salinity on optimized surface
%           ctns_i      conservative temperature on optimized surface
%           pns_i       pressure on optimized surface
%
% Calls:    mld.m, epsilon.m, var_on_surf.m
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%           gravitational acceleration  m/s^2
%           buoyancy frequency          s^-1
%           scale factors               m
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   Partially modified by P. Barker (2010-13)
%   Partially modified by S. Riha (2013)
%   Principal investigator: Trevor McDougall
%   type 'help analyze_surface' for more information
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

user_input;


% prepare data
%dbstop in mld at 50
cut_off_choice = mld(s,ct,p); % mixed-layer depth
   
breakout=false;
stop_wetting=false;
cnt_it_after_wetting=0;

% store indices of wetted grid points of the last three iterations. Necessary to avoid periodic
% cut-and-append behaviour with a period extending over multiple
% iterations (we have seen periods of 2 and assume here that 3 is the worst possible case)
iwetted_old={[],[],[]}; 
load('data/stationindex.mat');
% iterations of inversion
it=0; % start with it=0 and increment it after the initial surface is written to output
%while it<=nit_max;
while 1

    % diagnose
    if save_iterations;
        if it==0; % dummy values
            sx=nan; sy=nan; drho=nan; res=nan; b=nan; n2ns=nan; 
        end
        diagnose_and_write(it,sns,ctns,pns,sx,sy,drho,res,b,n2ns);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for comparison with gamma^i
    if breakout
        %rmsdrho=sqrt( nanmean( [drhodx(:);drhody(:)] .^2));
        rmsdrho=nan;
        
        % D_f
%         [nz,ny,nx]=size(s);
%         rho=gsw_rho(s(:,:),ct(:,:),p(:,:));
%         rho=reshape(rho,[nz,ny,nx]);
%         rhos=var_on_surf_stef(rho,p,pns);
%         
%         [n2,pmid]=gsw_Nsquared(s(:,:),ct(:,:),p(:,:));
%         n2=reshape(n2,[nz-1,ny,nx]);
%         pmid=reshape(pmid,[nz-1,ny,nx]);
%         n2s=var_on_surf_stef(n2,pmid,pns);
%         
%         fac=(1/9.81)*rhos.*n2s;
%         facx=0.5*(fac+circshift(fac,[0 -1])); % drhox staggered
%         facy=0.5*(fac+circshift(fac,[-1 0]));
%         
%         load('data/dxdy.mat') 
%         
%         sx=drhodx./(facx.*dx);
%         sy=drhody./(facy.*dy);
        s2_bar=nanmean( [sx(:);sy(:)] .^2);
        s2_med=nanmedian( [sx(:);sy(:)] .^2);
        mdf=1e3*s2_bar;
        df_med=1e3*s2_med;
        
        
%         if nx~=1
%             dx=0.5*(dx(:,1:end-1)+dx(:,2:end));
%             dx=horzcat(dx(:,end), dx); % sloppy 
%             ex=ex./dx;
%         end
%         dy=0.5*(dy(1:end-1,:)+dy(2:end,:));
%         dy=vertcat(dy(end,:), dy); % sloppy
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %if (it==nit_max) || breakout
    if breakout
        break
    end
    
    
    % Locations where outcropping occurs may have changed. Add points to
    % surface if necessary.
    %disp('Appending') 
    
    if ~stop_wetting
        [sns,ctns,pns,nneighbours,iwetted]=wetting(sns,ctns,pns,s,ct,p);
        %iwetted=[];
        %nneighbours=0;
    else
        nneighbours=0;
    end
    for ll=1:length(iwetted_old)
        iw_old=iwetted_old{ll};
        if size(iwetted)==size(iw_old)
            app_and_cut=all(iwetted==iw_old); % identical points keep getting appended and cut off
            if app_and_cut
                stop_wetting=true;
            end
        end
    end
    iwetted_old(1)=[]; % forget last
    iwetted_old{1,3}=iwetted; % insert most recent
    if nneighbours==0
        stop_wetting=true;
    end
    if it>5  
        if stop_wetting
            cnt_it_after_wetting=cnt_it_after_wetting+1;
        end
    end
    if cnt_it_after_wetting==nit_after_wetting
        breakout=true;
    end

    it=it+1; % start the next iteration
    
%    if (it==nit_max) && (nneighbours~=0)
%        error('need more wetting?')
%    end
    disp(['iter ',int2str(it), '    append ',num2str(nneighbours)]);
        
%     % disregard data above mixed layer depth
%     drhodx(pns<=cut_off_choice)=nan;
%     drhody(pns<=cut_off_choice)=nan;
%     sns(pns<=cut_off_choice)=nan;
%     ctns(pns<=cut_off_choice)=nan;
%     pns(pns<=cut_off_choice)=nan;

%    dbstop in n2 at 6

    if no_land_mask
        load('data/latlon.mat')
        [ocean, n] = gamma_ocean_and_n(s,ct,p,lon,lat);
        save('data/no_land_mask.mat', 'ocean', 'n')
    end
 
    % calculate delta^tilde rho
    %[drhodx,drhody]=delta_tilde_rho(sns,ctns,pns);
    [sx,sy]=slope_error(sns,ctns,pns,s,ct,p);
    %keyboard
    if use_b    
        [drhodx,drhody,regions,b]=use_bstar(drhodx,drhody,sns,ctns,pns,s,ct,p);
        %keyboard
    else
        %regions=find_regions(pns);
        [regions]=remove_points(sx,sy,pns);
%         setnan=true(size(sns));
%         regions=find_regions(sns);
%         for iregion=1:length(regions)
%             region=regions{iregion};
%             if ismember(istation,region)
%                 setnan(region)=false;
%                 pns(setnan)=nan;
%                 sns(setnan)=nan;
%                 ctns(setnan)=nan;
%             end
%         end
%         [sns,ctns,pns] = get_connected(sns,ctns,pns,istation);
    end
%keyboard
    %[drho,res]=solve_lsqr(regions, drhodx, drhody); 
    [drho,res]=solve_lsqr(regions, sx, sy); 
    setnan=true(size(sns));
    regions=find_regions(drho);
    for iregion=1:length(regions)
        region=regions{iregion};
        if ismember(istation,region)
            setnan(region)=false;
        end
    end
    drho(setnan)=nan;
    %keyboard
    if any(~isnan(drho(:)))

%         [nz,ny,nx]=size(s);
%         [n2,pmid]=gsw_Nsquared(s,ct,p);
%         n2=reshape(n2,[nz-1,ny,nx]);
%         pmid=reshape(pmid,[nz-1,ny,nx]);
%         n2s=var_on_surf_stef(n2,pmid,pns);
%         clip=1e-7;
%         n2s(n2s<clip)=clip;
%         rho=gsw_rho(s,ct,p);
%         rhos=var_on_surf_stef(rho,p,pns);
%         fac=9.81./(n2s.*rhos);
        
        %drho=drho./fac;
        disp(['max(drho): ',num2str(max(drho(:)))])
        pns=pns+0.1*drho;
        sns=var_on_surf_stef(s,p,pns);
        ctns=var_on_surf_stef(ct,p,pns);

        if use_b
            drho=drho./b;
        end

        %[sns, ctns, pns] = depth_ntp_simple(sns(:)', ctns(:)', pns(:)', s(:,:), ct(:,:), p(:,:), drho(:)' );

        [zi,yi,xi]=size(s);
        sns=reshape(sns,[yi xi]);
        ctns=reshape(ctns,[yi xi]);
        pns=reshape(pns,[yi xi]);

        % remove any regions which may have been detached during optimization
        %keyboard
        if clamp_on_point 
            load('data/stationindex.mat') % istation
            [sns,ctns,pns] = get_connected(sns,ctns,pns,istation);
        end
    else
        breakout=true;
    end

end

sns_i=sns;
ctns_i=ctns;
pns_i=pns;

end


%function [drhodx,drhody,sns,ctns,pns]=use_bstar(drhodx,drhody,sns,ctns,pns,s,ct,p)
function [drhodx,drhody,regions,b]=use_bstar(drhodx,drhody,sns,ctns,pns,s,ct,p)
    user_input;
    %keyboard
    [zi,yi,xi]=size(s);
    [n2,pmid]=gsw_Nsquared(s,ct,p);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    n2=reshape(n2,[zi-1,yi,xi]);
    pmid=reshape(pmid,[zi-1,yi,xi]);
   
    % regrid onto drhox and drhoy grids
    n2x=  0.5*(circshift(n2,   [0 0 -1])+n2);
    n2y=  0.5*(circshift(n2,   [0 -1 0])+n2);
    pmidx=0.5*(circshift(pmid, [0 0 -1])+pmid);
    pmidy=0.5*(circshift(pmid, [0 -1 0])+pmid);
 
    [rkx,rky]=grad_rhokappa(s,ct,p); % defined on drho grid
    rkx=0.5*(circshift(rkx,   [-1 0 0])+rkx); % regrid on vertical n2 grid
    rky=0.5*(circshift(rky,   [-1 0 0])+rky);
    rkx=rkx(1:end-1,:,:);
    rky=rky(1:end-1,:,:);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     n2=reshape(n2,[zi-1,yi,xi]);
%     n2=cat(1,n2,n2(end,:,:));
%     for ii=1:xi
%         for jj=1:yi
%             for kk=2:zi
%                 if isnan(n2(kk,jj,ii)) && ~isnan(s(kk,jj,ii))
%                     n2(kk,jj,ii)=n2(kk-1,jj,ii);
%                 end
%             end
%         end
%     end
%     
%     %keyboard
%     pmid=p;
%     %keyboard
%     % regrid onto drhox and drhoy grids
%     n2x=  0.5*(circshift(n2,   [0 0 -1])+n2);
%     n2y=  0.5*(circshift(n2,   [0 -1 0])+n2);
%     pmidx=0.5*(circshift(pmid, [0 0 -1])+pmid);
%     pmidy=0.5*(circshift(pmid, [0 -1 0])+pmid);
%  
%     [rkx,rky]=grad_rhokappa(s,ct,p); % defined on drho grid
% %     rkx=0.5*(circshift(rkx,   [-1 0 0])+rkx); % regrid on vertical n2 grid
% %     rky=0.5*(circshift(rky,   [-1 0 0])+rky);
% %     rkx=rkx(1:end-1,:,:);
% %     rky=rky(1:end-1,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    g=9.81;
    lnbx=(g^2./n2x).*rkx;
    lnby=(g^2./n2y).*rky;    
    
    pnsx= 0.5*(circshift(pns,    [0 -1])+pns);
    pnsy= 0.5*(circshift(pns,    [-1 0])+pns);
    
    lnbx=var_on_surf_stef(lnbx,pmidx,pnsx);
    lnby=var_on_surf_stef(lnby,pmidy,pnsy);

    
%     % regrid onto drhox and drhoy grids
%     n2x=  0.5*(circshift(n2,   [0 0 -1])+n2);
%     n2y=  0.5*(circshift(n2,   [0 -1 0])+n2);
%     pmidx=0.5*(circshift(pmid, [0 0 -1])+pmid);
%     pmidy=0.5*(circshift(pmid, [0 -1 0])+pmid);
%     pnsx= 0.5*(circshift(pns,    [0 -1])+pns);
%     pnsy= 0.5*(circshift(pns,    [-1 0])+pns);
%     
%     n2x_ns=var_on_surf_stef(n2x,pmidx,pnsx);
%     n2y_ns=var_on_surf_stef(n2y,pmidy,pnsy);
%     
%     [rkx,rky]=grad_rhokappa(s,ct,p); % defined on drho grid
%     
%     px=0.5*(circshift(p,[0 0 -1])+p);
%     py=0.5*(circshift(p,[0 -1 0])+p);
%     
%     rkx_ns=var_on_surf_stef(rkx,px,pnsx);
%     rky_ns=var_on_surf_stef(rky,py,pnsy);
%     
%     g=9.81;
%     lnbx=(g^2./n2x_ns.^2).*rkx_ns;
%     lnby=(g^2./n2y_ns.^2).*rky_ns;

    [regions]=remove_points(lnbx,lnby,pns);

    lnb=solve_lsqr(regions,lnbx,lnby);
    b=exp(lnb);
    
    bx=0.5*(circshift(b,[0 -1])+b);
    by=0.5*(circshift(b,[-1 0])+b);

    drhodx=bx.*drhodx;
    drhody=by.*drhody;
    
end

function [regions]=remove_points(lnbx,lnby,pns)
    goodx=~isnan(lnbx);
    goody=~isnan(lnby);
    
    % lnbx and lnby are staggered with respect to pns
    % tmp lives on the pns grid and is true if at least one (staggered) neighbour is not-nan.
    tmp=goodx | circshift(goodx,[0 1]);
    tmp= tmp | (goody | circshift(goody,[1 0]));
    
    % eqe is true if a tmp has an eastern neighbour...
    eqe=tmp & circshift(tmp,[0 -1]);
    eqn=tmp & circshift(tmp,[-1 0]);
    
    % ...there is an eastern neighbour, but not an eastward equation. Remove.
    remove_eq= eqe&~goodx;
    remove_eq=remove_eq | eqn&~goody;
    
    tmp(remove_eq)=false;
    
    tmp=double(tmp);
    tmp(tmp==0)=nan;
    
    disp(['      Lost ', num2str(sum(~isnan(pns(:)))-sum(~isnan(tmp(:)))), ...
              ' points during re-gridding'])
    if any(tmp(:))
        regions=find_regions(tmp);
    else
        regions={};
    end
end
    
function [sns,ctns,pns,nneighbours,iw]=wetting(sns,ctns,pns,s,ct,p)
% This function calls the actual wetting routine for each region
% separately, to avoid problems arising from two regions being separated by
% only one single wet point. In that case there could be wetting from only
% one of both directions (possibly the wrong one).
[yi,xi]=size(sns);

regions=find_regions(sns);
iw=false(xi*yi,1);
for ireg=1:length(regions)
    reg=regions{ireg};
    
    sns_r=nan*ones(yi,xi);
    ctns_r=nan*ones(yi,xi);
    pns_r=nan*ones(yi,xi);
    
    sns_r(reg)=sns(reg);
    ctns_r(reg)=ctns(reg);
    pns_r(reg)=pns(reg);
    

    [sns_r,ctns_r,pns_r,~,iw_r]=wetting_region(sns_r,ctns_r,pns_r,s,ct,p);

    sns(iw_r)=sns_r(iw_r);
    ctns(iw_r)=ctns_r(iw_r);
    pns(iw_r)=pns_r(iw_r);
    iw=iw|iw_r;
    
end

nneighbours=sum(iw);

end

function [sns,ctns,pns,nneighbours,iwetted]=wetting_region(sns,ctns,pns,s,ct,p)

user_input;

[yi,xi]=size(sns);

wet=~isnan(squeeze(s(1,:,:))) & isnan(sns); % wet points at ocean surface excluding ans
wets=~isnan(sns); % wet points on ans
nn=wet(:) & circshift(wets(:),-1); % wet points with north. neighbour on ans
sn=wet(:) & circshift(wets(:),1);
en=wet(:) & circshift(wets(:),-yi);
wn=wet(:) & circshift(wets(:),yi);

nn(yi:yi:yi*xi)=false;
sn(1:yi:(xi-1)*yi+1)=false;
if ~zonally_periodic;
    en((xi-1)*yi+1:xi*yi)=false;
    wn(1:yi)=false;
end

% if a point adjacent to ans boundary has multiple neighbours, just do one neutral
% calculation
% TODO: start in eastward direction? Trevor has preference, see notes.
wn=wn & ~en;
nn=nn & ~wn & ~en;
sn=sn & ~nn & ~wn & ~en;

inds=[1:xi*yi]';
inds_neighbour=circshift(inds,-yi);
neighbour=inds_neighbour(en);
[sns(en),ctns(en),pns(en)] = depth_ntp_simple(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,en),ct(:,en),p(:,en)); 


inds_neighbour=circshift(inds,yi);
neighbour=inds_neighbour(wn);
[sns(wn),ctns(wn),pns(wn)] = depth_ntp_simple(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,wn),ct(:,wn),p(:,wn)); 

inds_neighbour=circshift(inds,-1);
neighbour=inds_neighbour(nn);
[sns(nn),ctns(nn),pns(nn)] = depth_ntp_simple(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,nn),ct(:,nn),p(:,nn)); 

inds_neighbour=circshift(inds,1);
neighbour=inds_neighbour(sn);
[sns(sn),ctns(sn),pns(sn)] = depth_ntp_simple(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,sn),ct(:,sn),p(:,sn)); 

s1=sum(~isnan(sns(en)));
s2=sum(~isnan(sns(wn)));
s3=sum(~isnan(sns(nn)));
s4=sum(~isnan(sns(sn)));

nneighbours=s1+s2+s3+s4;

iwetted= ~isnan(sns(:)) & (en | wn | sn | nn);
if sum(iwetted)~=nneighbours
    error('something is wrong')
end
%disp(['Number of points added: ',num2str(nneighbours)])

end


function [drho,res]=solve_lsqr(regions, xx, yy)
user_input;
if no_land_mask
    load('data/no_land_mask.mat')
end
if clamp_on_point
    load('data/stationindex.mat') % istation
end
[yi,xi]=size(xx);
drho = nan(yi,xi);

for iregion=1:length(regions)
    
    region=regions{iregion};
    
    % set up east-west equations for weighted inversion
    reg=false(1,xi*yi)';
    reg(region)=true;
    
    en= reg & circshift(reg,-yi); %  find points between which a zonal gradient can be computed. en is true at a point if its eastward neighbor is in the region

    if no_land_mask % in case there is no land mask (e.g. in a climatology data set with too coarse resolution to represent land)
        % 'en_in_other_basin' is true if there is land (in reality) between this point and its eastern neighbour
        error('find_regions() doesn''t have this functionality yet')
        test=ocean(:).*circshift(ocean(:),-yi);
        en_in_other_basin= (test==5) & reg;
        en(en_in_other_basin)=false;
    end
    if ~zonally_periodic;  % remove equations for eastern boundary for zonally-nonperiodic domain
        en((xi-1)*yi+1:xi*yi)=false;
    end
    sreg=cumsum(reg); % sparse indices of points forming the region (points of non-region are indexed with dummy)
    sreg_en=circshift(sreg,-yi); % sparse indices of eastward neighbours
    
    j1_ew=sreg_en(en);  % j1 are j-indices for matrix coefficient 1
    j2_ew=sreg(en); % j2 are j-indices for matrix coefficient -1
    
    % set up north-south equations for weighted inversion
    nn= reg & circshift(reg,-1);
    if no_land_mask % in case there is no land mask (e.g. in a climatology data set with too coarse resolution to represent land)
        % 'nn_in_other_basin' is true if there is land (in reality) between this point and its northern neighbour
        test=ocean(:).*circshift(ocean(:),-1);
        nn_in_other_basin= (test==5) & reg;
        nn(nn_in_other_basin)=false;
    end
    nn(yi:yi:xi*yi)=false; % remove equations for northern boundary
    sreg_nn=circshift(sreg,-1);
    
    j1_ns=sreg_nn(nn);
    j2_ns=sreg(nn);
    
    % clamp at point
    if ismember(istation,region)
        j1_condition=sreg(istation);
    else
        j1_condition=[1:sum(reg)];
    end
    
    j1=[j1_ew',j1_ns',j1_condition];
    j2=[j2_ew',j2_ns'];
    
    i2=1:(sum(en)+sum(nn)); % i-indices for matrix coeff. -1
    if ismember(istation,region)
        i1=[i2, (sum(en)+sum(nn)+1)]; % i-indices for matrix coeff. 1
    else
        i1=[i2, (sum(en)+sum(nn)+1)*ones(1,sum(reg))]; % i-indices for matrix coeff. 1
    end

    
    % build sparse matrices
    A=sparse([i1,i2],[j1,j2],[ones(1,length(i1)),-ones(1,length(i2))]);
    b=sparse( [xx(en); yy(nn); 0 ]);
    
    %disp(['solving for region ',int2str(iregion)]);
    switch solver
        case 'iterative'
            tic
            [x,flag,relres,iter,resvec,lsvec] = lsqr(A,b,1e-15,10000);
            %display(['LSQR() took ',num2str(toc),' seconds for ',num2str(length(lsvec)),' iterations']);
            res=relres;
%             if flag ~=0 
%                 disp('LSQR did not converge')
%                 if flag~=1
%                     disp('MAXIT WAS NOT REACHED')
%                 end
%             end
        case 'exact'
            x = (A'*A)\(A'*b);
    end
    
    x = full(x)';
    
    % put density changes calculated by the least-squares solver into
    % their appropriate position in the matrix

    drho(region) = x;
    
end
end



function [sns_out,ctns_out,pns_out] = dz_from_drho(sns, ctns, pns, s, ct, p, drho );

user_input; % read delta
[zi,yi,xi]=size(s);

rho_surf=gsw_rho(sns(:),ctns(:),pns(:));
t2=rho_surf-drho(:);

inds=1:yi*xi;
fr=true(1,yi*xi);
pns_out = nan(yi,xi);
sns_out = nan(yi,xi);
ctns_out = nan(yi,xi);
pns=pns(:);

% % don't adjust at locations where abs(drho) is almost zero
% noadjust=abs(drho(:))<=delta;
% fr(noadjust)=false;
% inds=inds(fr);
% sns_out(noadjust)=sns(noadjust);
% ctns_out(noadjust)=ctns(noadjust);
% pns_out(noadjust)=pns(noadjust);
% s=s(:,fr);
% ct=ct(:,fr);
% p=p(:,fr);

pns_stacked=repmat(pns(fr)',[zi 1]); % stack pressure of current surface vertically
t2_stacked=repmat(t2(fr)',[zi 1]); % stack locally referenced density of current surface vertically

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative code

% stack=zi;
% ii=bsxfun(@times,1:yi*xi,ones(stack,1));
% pns_stacked=pns(ii); % stack pressure of current surface vertically
% t2_stacked_full=t2(ii); % stack locally referenced density of current surface vertically
% 
% t1=gsw_rho(s(:,:),ct(:,:),pns_stacked(:,inds)); % 3-d density referenced to pressure of the current surface
% t2_3d=t2_stacked_full(:,inds); %
% tni=t1-t2_3d; % rho-(rho_s+rho'); find corrected surface by finding roots of this term
% 
% Itni_p = tni>0;
% Itni_n = tni<0;
% 
% zc = any(Itni_p,1) & any(Itni_n,1); % horizontal indices of locations where zero-crossing occurs
% [min_tni, lminr]=min(abs(tni));
% cond1=min_tni>delta;
% final=min_tni<=delta; % at these horizontal indices root has been found
% fr= zc & cond1; %  at these horizontal locations we have to use root finding
% 
% lminr=lminr+stack*[0:size(Itni_n,2)-1];
% lminr=lminr(final);
% 
% sns_out(inds(final))=s(lminr); % adjust surface where root has already been found
% ctns_out(inds(final)) =ct(lminr);
% pns_out(inds(final)) =p(lminr);
% 
% k=sum(Itni_n,1); % find indices of flattened 3d-array where we use root finding
% k=k+stack*[0:size(Itni_n,2)-1];
%
% options=optimset('TolX',1e-3);
% for ii=inds(fr);
%     su=s(k(ii)); sl=s(k(ii)+1); % linear interpolation of s and ct
%     ctu=ct(k(ii)); ctl=ct(k(ii)+1);
%     pu=p(k(ii)); pl=p(k(ii)+1);
%     
%     sp=@(p) su+(sl-su)*(p-pu)/(pl-pu);
%     ctp=@(p) ctu+(ctl-ctu)*(p-pu)/(pl-pu);
%     fac=1./sqrt(gsw_rho(su,ctu,pu)^2+gsw_rho(sl,ctl,pl)^2); % scale to avoid having to set tolerance ?
%     drho_normalized=@(p)  fac.*(gsw_rho(sp(p), ctp(p), pns(ii)) -t2(ii));
%     
%     %dmyrho_normalized=@(p,su,sl,ctu,ctl,pu,pl)  1.*(gsw_rho( su+(sl-su)*(p-pu)/(pl-pu), ctu+(ctl-ctu)*(p-pu)/(pl-pu), pns(ii))-t2(ii));
%     %drho_normalized=@(p) dmyrho_normalized(p,su,sl,ctu,ctl,pu,pl);
%     proot=fzero(drho_normalized, [pu,pl],options);
%     
%     sns_out(ii)= su+(sl-su)*(proot-pu)/(pl-pu);
%     ctns_out(ii)=ctu+(ctl-ctu)*(proot-pu)/(pl-pu);
%     pns_out(ii)=proot;
% end

% alternative code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative code
    
refine_ints=100;

cnt=0;
while 1
    cnt=cnt+1;
    
    t1=gsw_rho(s(:,:),ct(:,:),pns_stacked); % 3-d density referenced to pressure of the current surface

    F=t1-t2_stacked; % rho-(rho_s+rho'); find corrected surface by finding roots of this term
    
    %dbstop in root_core at 11
    [s,ct,p,sns_out,ctns_out,pns_out, inds, fr]=root_core(F,inds,refine_ints,s,ct,p,sns_out,ctns_out,pns_out);
    
    if all(~fr) % break out of loop if all roots have been found
        break
    end
    
    pns_stacked=pns_stacked(1:refine_ints+1,fr);
    t2_stacked=t2_stacked(1:refine_ints+1,fr);
    
end

% alternative code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



function diagnose_and_write(it,sns,ctns,pns,drhodx,drhody,drho,res,b,n2ns)
user_input; % read nit, etc.
load('data/dxdy.mat') % for epsilon
if it==0 % initialize
    [yi,xi]=size(sns);
    
    sns_hist = nan(1,yi,xi); 
    ctns_hist = nan(1,yi,xi);
    pns_hist = nan(1,yi,xi);
    
    %slope_square = nan(nit,1);
    %eps_rms_hist=nan(nit,1);
    
    drhoxy_rms_hist=nan(1,1); 
    drho_rms_hist=nan(1,1);
    epsilon_rms_hist=nan(1,1);
    
    drho_hist = nan(1,yi,xi);
    res_hist = nan(1,1);
    
    b_hist = nan(1,yi,xi);
    n2ns_hist = nan(1,yi,xi);
    
    vars = {'sns_hist','ctns_hist','pns_hist','drhoxy_rms_hist','drho_rms_hist','epsilon_rms_hist','drho_hist','res_hist','b_hist','n2ns_hist'};
    save(history_file, vars{:},'-v7.3');
end

iteration_history = matfile(history_file,'Writable',true);
iteration_history.sns_hist(it+1,:,:) = permute(sns,[3 1 2]);
iteration_history.ctns_hist(it+1,:,:) = permute(ctns,[3 1 2]);
iteration_history.pns_hist(it+1,:,:) = permute(pns,[3 1 2]);

if it>0
    iteration_history.drho_hist(it,:,:) = permute(drho,[3,1,2]);
    iteration_history.res_hist(it,1) = res;
    
    iteration_history.b_hist(it,:,:) = permute(b,[3,1,2]);
    iteration_history.n2ns_hist(it,:,:) = permute(n2ns,[3,1,2]);
    
    %s1=ex(~isnan(ex));  
    %s2=ey(~isnan(ey));
    %square=[s1(:) ; s2(:)].^2;
    %slope_square(it,1) = sum(square);
    %no_pts =length(square);
    %iteration_history.eps_rms_hist(it,1) = sqrt(slope_square(it,1)/no_pts);
    
    iteration_history.drho_rms_hist(it,1)= sqrt( nanmean(drho(:).^2) );
    iteration_history.drhoxy_rms_hist(it,1)= sqrt( nanmean( [drhodx(:);drhody(:)] .^2)); % staggerd grid
    iteration_history.epsilon_rms_hist(it,1)= sqrt( nanmean( [drhodx(:)./dx(:);drhody(:)./dy(:)] .^2)); % staggerd grid
end

end





