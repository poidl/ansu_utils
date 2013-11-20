function [sns_i,ctns_i,pns_i] = optimize_surface_exact(s,ct,p,sns,ctns,pns)

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

stop_wetting=false;

% iterations of inversion
it=0; % start with it=0 and increment it after the initial surface is written to output
while it<=nit;
    
    % diagnose
    if save_iterations;
        if it==0; % dummy values
            drhodx=nan; drhody=nan; drho=nan; b=nan; n2ns=nan;
        end
        diagnose_and_write(it,sns,ctns,pns,drhodx,drhody,drho,b,n2ns);
    end
    
    if it==nit; % break out if done
        break
    end
    
    it=it+1; % start the next iteration
    disp(['Iteration nr.',int2str(it)]);
    

%     % Locations where outcropping occurs may have changed. Add points to
%     % surface if necessary.
%     if it<nit && ~stop_wetting;
%         disp('Wetting')
%         if (it==1||it==2); % it==2 may not be necessary for good starting surfaces, but it is necessary when starting from an isobar
%         %if (it==1); 
%             [sns,ctns,pns,nneighbours]=wetting(sns,ctns,pns,s,ct,p);
%         else
%             nneighbours_old=nneighbours;
%             [sns,ctns,pns,nneighbours]=wetting(sns,ctns,pns,s,ct,p);
%             if nneighbours>=nneighbours_old;
%                 stop_wetting=true;
%                 disp(['Number of wet points added is equal or larger to previous iteration. Stop wetting.'])
%             end
%         end
%     end

    
    % calculate delta^tilde rho
    [drhodx,drhody]=delta_tilde_rho(sns,ctns,pns);
        
    
%     % disregard data above mixed layer depth
%     drhodx(pns<=cut_off_choice)=nan;
%     drhody(pns<=cut_off_choice)=nan;
%     sns(pns<=cut_off_choice)=nan;
%     ctns(pns<=cut_off_choice)=nan;
%     pns(pns<=cut_off_choice)=nan;

%    dbstop in n2 at 6
    if 1
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

        ro=gsw_rho(sns,ctns,pns);
        rox=0.5*(circshift(ro, [0 -1])+ro);
        if ~zonally_periodic;
            rox(:,xi) = nan;
        end
        roy=0.5*(circshift(ro, [-1 0])+ro);
        roy(yi,:) = nan;

        [kx_ns,ky_ns]=kappaxy(s,ct,p,pns);

        fx=9.81^2*rox.*(1./n2nsx).*kx_ns;
        fy=9.81^2*roy.*(1./n2nsy).*ky_ns;

        iyeq=(~isnan(pns) & ~isnan(circshift(pns,[-1 0])));
        bad=(iyeq & isnan(ky_ns));

        pns(bad)=nan;
        sns(bad)=nan;
        ctns(bad)=nan;

        % find independent regions -> a least-squares problem is solved for
        % each of these regions
        regions=find_regions(pns);

        lnb=solve_lsqr(regions,fx,fy);
        b=exp(lnb);

        bx=0.5*(circshift(b,[0 -1])+b);
        by=0.5*(circshift(b,[-1 0])+b);
        
        drhodx_mod=bx.*drhodx;
        drhody_mod=by.*drhody;
        regions=find_regions(pns);
 
        % solve for delta rho
        drho=solve_lsqr(regions, drhodx_mod, drhody_mod);  
        drho=drho./b;

    else
    
        regions=find_regions(pns);

        % solve for delta rho
        drho=solve_lsqr(regions, drhodx, drhody);    
    
    end
    
    % find corrected surface
    [sns, ctns, pns] = dz_from_drho(sns, ctns, pns, s, ct, p, drho );
    
end

sns_i=sns;
ctns_i=ctns;
pns_i=pns;

end


function [sns,ctns,pns,nneighbours]=wetting(sns,ctns,pns,s,ct,p)

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

[sns(en),ctns(en),pns(en)] = depth_ntp_iter(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,en),ct(:,en),p(:,en)); 

inds_neighbour=circshift(inds,yi);
neighbour=inds_neighbour(wn);
[sns(wn),ctns(wn),pns(wn)] = depth_ntp_iter(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,wn),ct(:,wn),p(:,wn)); 

inds_neighbour=circshift(inds,-1);
neighbour=inds_neighbour(nn);
[sns(nn),ctns(nn),pns(nn)] = depth_ntp_iter(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,nn),ct(:,nn),p(:,nn)); 

inds_neighbour=circshift(inds,1);
neighbour=inds_neighbour(sn);
[sns(sn),ctns(sn),pns(sn)] = depth_ntp_iter(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,sn),ct(:,sn),p(:,sn)); 

s1=sum(~isnan(sns(en)));
s2=sum(~isnan(sns(wn)));
s3=sum(~isnan(sns(nn)));
s4=sum(~isnan(sns(sn)));

nneighbours=s1+s2+s3+s4;
disp(['Number of points added: ',num2str(nneighbours)])

end




function drho=solve_lsqr(regions, xx, yy)
user_input;
[yi,xi]=size(xx);
drho = nan(yi,xi);

for iregion=1:length(regions)
    
    region=regions{iregion};
    
    % set up east-west equations for weighted inversion
    reg=false(1,xi*yi)';
    reg(region)=true;
    
    en= reg & circshift(reg,-yi); %  find points between which a zonal gradient can be computed. en is true at a point if its eastward neighbor is in the region
    if ~zonally_periodic;  % remove equations for eastern boundary for zonally-nonperiodic domain
        en((xi-1)*yi+1:xi*yi)=false;
    end
    sreg=cumsum(reg); % sparse indices of points forming the region (points of non-region are indexed with dummy)
    sreg_en=circshift(sreg,-yi); % sparse indices of eastward neighbours
    
    j1_ew=sreg_en(en);  % j1 are j-indices for matrix coefficient 1
    j2_ew=sreg(en); % j2 are j-indices for matrix coefficient -1
    
    % set up north-south equations for weighted inversion
    nn= reg & circshift(reg,-1);
    nn(yi:yi:xi*yi)=false; % remove equations for northern boundary
    sreg_nn=circshift(sreg,-1);
    
    j1_ns=sreg_nn(nn);
    j2_ns=sreg(nn);
    
    % make the average of the potential zero
    % this should keep the surface from drifting away from the initial condition
    % we might change that to a different condition
    j1_condition=[1:sum(reg)];
    
    j1=[j1_ew',j1_ns',j1_condition];
    j2=[j2_ew',j2_ns'];
    
    i2=1:(sum(en)+sum(nn)); % i-indices for matrix coeff. -1
    i1=[i2, (sum(en)+sum(nn)+1)*ones(1,sum(reg))]; % i-indices for matrix coeff. 1
    
    % build sparse matrices
    A=sparse([i1,i2],[j1,j2],[ones(1,length(i1)),-ones(1,length(i2))]);
    b=sparse( [xx(en); yy(nn); 0 ]);
    
    disp(['solving for region ',int2str(iregion)]);
    switch solver
        case 'iterative'
            [x,dummyflag] = lsqr(A,b,delta,1000);
        case 'exact'
            x = (A'*A)\(A'*b);
    end
    
    x = full(x)';
    
    % put density changes calculated by the least-squares solver into
    % their appropriate position in the matrix

    drho(region) = x;
    
end
end







function diagnose_and_write(it,sns,ctns,pns,drhodx,drhody,drho,b,n2ns)
user_input; % read nit, etc.

if it==0 % initialize
    [yi,xi]=size(sns);
    
    sns_hist = nan(nit+1,yi,xi); % store variables on initial surface and on nit improvements (=> nit+1)
    ctns_hist = nan(nit+1,yi,xi);
    pns_hist = nan(nit+1,yi,xi);
    
    %slope_square = nan(nit,1);
    %eps_rms_hist=nan(nit,1);
    
    drhoxy_rms_hist=nan(nit,1); 
    drho_rms_hist=nan(nit,1);
    
    drho_hist = nan(nit,yi,xi);
    
    b_hist = nan(nit,yi,xi);
    n2ns_hist = nan(nit,yi,xi);
    
    vars = {'sns_hist','ctns_hist','pns_hist','drhoxy_rms_hist','drho_rms_hist','drho_hist','b_hist','n2ns_hist'};
    save(history_file, vars{:},'-v7.3');
end

iteration_history = matfile(history_file,'Writable',true);
iteration_history.sns_hist(it+1,:,:) = permute(sns,[3 1 2]);
iteration_history.ctns_hist(it+1,:,:) = permute(ctns,[3 1 2]);
iteration_history.pns_hist(it+1,:,:) = permute(pns,[3 1 2]);

if it>0
    iteration_history.drho_hist(it,:,:) = permute(drho,[3,1,2]);
    
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
    
end

end





