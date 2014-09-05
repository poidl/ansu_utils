function [sns_i,ctns_i,pns_i,errx,erry] = optimize_surface_exact(s,ct,p,sns,ctns,pns)

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

%cut_off_choice = mld(s,ct,p); % mixed-layer depth
   
breakout=false;
stop_wetting=false;
cnt_it_after_wetting=0;

% store indices of wetted grid points of the last three iterations. Necessary to avoid periodic
% cut-and-append behaviour with a period extending over multiple
% iterations (we have seen periods of 2 and assume here that 3 is the worst possible case)
iwetted_old={[],[],[]}; 

% iterations of inversion
it=0; % start with it=0 and increment it after the initial surface is written to output
%while it<=nit_max;
while 1

    % diagnose
    if save_iterations;
        if it==0; % dummy values
            errx=nan; erry=nan; derr=nan; res=nan; b=nan; n2ns=nan; 
        end
        diagnose_and_write(it,sns,ctns,pns,errx,erry,derr,res,b,n2ns);
    end
    
    %if (it==nit_max) || breakout
    if breakout
        break
    end
    
    % Locations where outcropping occurs may have changed. Add points to
    % surface if necessary.
    %disp('Appending') 
    
    if ~stop_wetting
        [sns,ctns,pns,nneighbours,iwetted]=wetting(sns,ctns,pns,s,ct,p);
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

    if no_land_mask
        load('data/latlon.mat')
        [ocean, n] = gamma_ocean_and_n(s,ct,p,lon,lat);
        save('data/no_land_mask.mat', 'ocean', 'n')
    end
 
    if strcmp(error_measure,'drho_local')
        [errx,erry]=delta_tilde_rho(sns,ctns,pns);
    elseif strcmp(error_measure,'slope_difference')
        [errx,erry]=slope_error(sns,ctns,pns,s,ct,p);
    end
    
    if use_b    
        if strcmp(error_measure,'drho_local')
            [errx,erry,regions,b]=use_bstar(errx,erry,pns,s,ct,p);
            [derr,res]=solve_lsqr(regions, errx, erry);
        else 
            error('not implemented')
        end
    else
        if strcmp(error_measure,'drho_local')
            regions=find_regions(pns);
            [derr,res]=solve_lsqr(regions, errx, erry); 
            
        elseif strcmp(error_measure,'slope_difference')
            regions=find_regions(pns);
%             [regions]=remove_points(errx,erry,pns);
             [derr,res]=solve_lsqr(regions, errx, erry);
%             % only keep region at backbone
%             setnan=true(size(sns));
%             regions=find_regions(derr);
%             for iregion=1:length(regions)
%                 region=regions{iregion};
%                 if ismember(istation,region)
%                     setnan(region)=false;
%                 end
%             end
%             derr(setnan)=nan;
        end
    end 
    
    if use_b && strcmp(error_measure,'drho_local')
        derr=derr./b;
    end
    
    if strcmp(error_measure,'drho_local')
        [sns, ctns, pns] = depth_ntp_simple(sns(:)', ctns(:)', pns(:)', s(:,:), ct(:,:), p(:,:), derr(:)' );
        [zi,yi,xi]=size(s);
        sns=reshape(sns,[yi xi]);
        ctns=reshape(ctns,[yi xi]);
        pns=reshape(pns,[yi xi]);
    elseif strcmp(error_measure,'slope_difference')
        pns=pns+0.1*derr;
        sns=var_on_surf_stef(s,p,pns);
        ctns=var_on_surf_stef(ct,p,pns);
    end
       
    % remove any regions which may have been detached during optimization
    %keyboard
    if clamp_on_point
        load('data/stationindex.mat') % istation
        [sns,ctns,pns] = get_connected(sns,ctns,pns,istation);
    end

end

sns_i=sns;
ctns_i=ctns;
pns_i=pns;

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


function [drhodx,drhody,regions,b]=use_bstar(drhodx,drhody,pns,s,ct,p)
    user_input;
    %keyboard
    
    g=9.81;
    
    [zi,yi,xi]=size(s);
    [n2,pmid]=gsw_Nsquared(s,ct,p);
        
    n2=reshape(n2,[zi-1,yi,xi]);
    pmid=reshape(pmid,[zi-1,yi,xi]);
 
    % regrid onto drhox and drhoy grids
    n2x=  0.5*(circshift(n2,   [0 0 -1])+n2);
    n2y=  0.5*(circshift(n2,   [0 -1 0])+n2);
    pmidx=0.5*(circshift(pmid, [0 0 -1])+pmid);
    pmidy=0.5*(circshift(pmid, [0 -1 0])+pmid);
 
    [rkx,rky]=delta_rhokappa(s,ct,p); 
    
    rkx=0.5*(circshift(rkx,   [-1 0 0])+rkx); % regrid onto vertical n2 grid
    rky=0.5*(circshift(rky,   [-1 0 0])+rky);
    rkx=rkx(1:end-1,:,:);
    rky=rky(1:end-1,:,:);

    lnbx=(g^2./n2x).*rkx;
    lnby=(g^2./n2y).*rky;    
    
    pnsx= 0.5*(circshift(pns,    [0 -1])+pns);
    pnsy= 0.5*(circshift(pns,    [-1 0])+pns);
    
    lnbx=var_on_surf_stef(lnbx,pmidx,pnsx);
    lnby=var_on_surf_stef(lnby,pmidy,pnsy);

    [regions]=remove_points(lnbx,lnby,pns);

    lnb=solve_lsqr(regions,lnbx,lnby);
    b=exp(lnb);
    
    bx=0.5*(circshift(b,[0 -1])+b);
    by=0.5*(circshift(b,[-1 0])+b);

    drhodx=bx.*drhodx;
    drhody=by.*drhody;
    
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


function [derr,res]=solve_lsqr(regions, xx, yy)
user_input;
if no_land_mask
    load('data/no_land_mask.mat')
end
if clamp_on_point
    load('data/stationindex.mat') % istation
end
[yi,xi]=size(xx);
derr = nan(yi,xi);

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

    derr(region) = x;
    
end
end


function diagnose_and_write(it,sns,ctns,pns,errx,erry,derr,res,b,n2ns)
user_input; % read nit, etc.
load('data/dxdy.mat') % for epsilon
if it==0 % initialize
    [yi,xi]=size(sns);
    
    sns_hist = nan(1,yi,xi); 
    ctns_hist = nan(1,yi,xi);
    pns_hist = nan(1,yi,xi);
    
    %slope_square = nan(nit,1);
    %eps_rms_hist=nan(nit,1);
    
    derrxy_rms_hist=nan(1,1); 
    derr_rms_hist=nan(1,1);
    epsilon_rms_hist=nan(1,1);
    
    derr_hist = nan(1,yi,xi);
    res_hist = nan(1,1);
    
    b_hist = nan(1,yi,xi);
    n2ns_hist = nan(1,yi,xi);
    
    vars = {'sns_hist','ctns_hist','pns_hist','derrxy_rms_hist','derr_rms_hist','epsilon_rms_hist','derr_hist','res_hist','b_hist','n2ns_hist'};
    save(history_file, vars{:},'-v7.3');
end

iteration_history = matfile(history_file,'Writable',true);
iteration_history.sns_hist(it+1,:,:) = permute(sns,[3 1 2]);
iteration_history.ctns_hist(it+1,:,:) = permute(ctns,[3 1 2]);
iteration_history.pns_hist(it+1,:,:) = permute(pns,[3 1 2]);

if it>0
    iteration_history.derr_hist(it,:,:) = permute(derr,[3,1,2]);
    iteration_history.res_hist(it,1) = res;
    
    iteration_history.b_hist(it,:,:) = permute(b,[3,1,2]);
    iteration_history.n2ns_hist(it,:,:) = permute(n2ns,[3,1,2]);
    
    %s1=ex(~isnan(ex));  
    %s2=ey(~isnan(ey));
    %square=[s1(:) ; s2(:)].^2;
    %slope_square(it,1) = sum(square);
    %no_pts =length(square);
    %iteration_history.eps_rms_hist(it,1) = sqrt(slope_square(it,1)/no_pts);
    
    iteration_history.derr_rms_hist(it,1)= sqrt( nanmean(derr(:).^2) );
    iteration_history.derrxy_rms_hist(it,1)= sqrt( nanmean( [errx(:);erry(:)] .^2)); % staggerd grid
    iteration_history.epsilon_rms_hist(it,1)= sqrt( nanmean( [errx(:)./dx(:);erry(:)./dy(:)] .^2)); % staggerd grid
end

end





