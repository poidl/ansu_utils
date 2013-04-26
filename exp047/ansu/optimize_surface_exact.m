function [sns_i,ctns_i,pns_i,ithist] = optimize_surface_exact(s,ct,p,g,n2,sns,ctns,pns,e1t,e2t,settings)

%           Optimize density surfaces to minimise the fictitious diapycnal diffusivity
%
% Usage:    [sns_i,ctns_i,pns_i] =
%           Optimize_surface(s,ct,p,g,n2,sns,ctns,pns,e1t,e2t,nit,choice,wrap)
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
%           g           gravitational acceleration
%           n2          buoyancy frequency
%           sns         salinity on initial density surface
%           ctns        conservative temperature on initial density
%                       surface
%           pns         pressure on initial density surface
%           e1t         grid length in meters in x-direction
%           e2t         grid length in meters in y-direction
%
% %           settings    set in main window.
% %           nit         maximum number of iterations is set to 150
% %           choice      choice between
% %                       's' slope error
% %                       'epsilon' density gradient error
% %           wrap        'none'
% %                       'long'
% %           solver      'iterative'
% %                       'exact'
% Output:   sns_i       salinity on optimized surface
%           ctns_i      conservative temperature on optimized surface
%           pns_i       pressure on optimized surface
%
% Calls:    mld.m, slope_error.m, var_on_surf.m
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
%   type 'help analyze_surface' for more information
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

wrap = settings.wrap;
choice = settings.slope;
solver = settings.solver;
nit = settings.nit;



[zi,yi,xi] = size(s);

%% calculate buoyancy frequency on density surface

p_mid = 0.5*(p(2:zi,:,:) + p(1:zi-1,:,:));


%% prepare data
slope_square = nan(nit+1,1);

sns_hist = nan(nit+1,yi,xi); % store variables on initial surface and on nit improvements (=> nit+1)
ctns_hist = nan(nit+1,yi,xi);
pns_hist = nan(nit+1,yi,xi);
eps_ss=nan(nit+1,1);

phiprime_e_hist = nan(nit,yi,xi);

cut_off_choice(1,1:yi,1:xi) = mld(s,ct,p);


%% iterations of inversion
    
it=0; % start with it=0 and increment it after the initial surface is written to output
while it<=nit;
    
    n2_ns = var_on_surf(pns,p_mid,n2);

    %% calculate slope errors/density gradient errors on initial density surface

    [ss,sx,sy,curl_s,ee,ex,ey,curl_e,ver] = slope_error(p,g,n2,sns,ctns,pns,e1t,e2t,'bp',wrap); %#ok
    
    
    %% diagnose
    switch choice
        case 'epsilon'
            square = ex .* ex + ey .* ey;
            slope_square(it+1,1) = nansum(square(:));
            no_pts = length(find(~isnan(square(:))));
            eps_ss(it+1,1) = sqrt(slope_square(it+1,1)/no_pts);
            
        case 's'
            square = sx .* sx + sy .* sy;
            slope_square(it+1,1) = nansum(square(:));
    end
   
    % store history
     
    sns_hist(it+1,:,:) = squeeze(sns);
    ctns_hist(it+1,:,:) = squeeze(ctns);
    pns_hist(it+1,:,:) = squeeze(pns(1,:,:));
    
    if it>0;
        phiprime_e_hist(it,:,:) = squeeze(phiprime_e);
    end
    
    if it==nit;
        break
    end
    
    it=it+1;
    disp(['Iteration nr.',int2str(it)]); % print number of iteration

    %% choice between minimizing slope errors or density gradient errors

    switch choice
        case 's'
            xx = sx;
            yy = sy;
        case 'epsilon'
            xx = ex;
            yy = ey;
    end
    
    if it>0;
        xx(pns<= cut_off_choice)=nan;
        yy(pns<= cut_off_choice)=nan;
    end
    
    xx_squeeze = squeeze(xx);
    yy_squeeze = squeeze(yy);

    %% disregard data above mixed layer depth
    
    n2_ns(pns<=cut_off_choice)=nan;
    pns(pns<=cut_off_choice)=nan;
    
    pns = squeeze(pns);
    n2_ns = squeeze(n2_ns);


    phiprime = nan(yi,xi);
    phiprime_e = nan(yi,xi);
    
    % flag wet points
    
    wet=~isnan(n2_ns);
    
    % find independent regions -> a least-squares problem is solved for
    % each of these regions
    
    cc=bwconncomp(wet,4);
    
    if strcmp(wrap,'long') % in a zonally periodic domain merge regions which are cut apart by grid boundary

        % find regions which have points at western and eastern boundary
        bdy_wet=false(1,length(cc.PixelIdxList));
        ireg=1:length(cc.PixelIdxList);
        for ii=ireg
            if any(cc.PixelIdxList{ii}<=yi | cc.PixelIdxList{ii}>yi*(xi-1))
                bdy_wet(ii)=true;
            end
        end
        iw=ireg(bdy_wet);
        
        merged=false(1,length(cc.PixelIdxList));

        ii=1;
        while ii<=length(iw)
            for jj=1: length(iw)
                if ii~=jj && ~(merged(iw(ii)) || merged(iw(jj)))
                    pts1=cc.PixelIdxList{iw(ii)};
                    pts2=cc.PixelIdxList{iw(jj)};
                    % check if western border of region iw(ii) intersects
                    % with eastern border of region iw(jj), or vice versa
                    cond1=~isempty( intersect( pts1(pts1<=yi)+yi*(xi-1),pts2(pts2>yi*(xi-1)) ) );
                    cond2=~isempty( intersect( pts2(pts2<=yi)+yi*(xi-1),pts1(pts1>yi*(xi-1)) ) );
                    if cond1 | cond2
                        cc.PixelIdxList{iw(ii)}=union( pts1, pts2 );
                        merged(iw(jj))=true; % iw(jj) has been merged into iw(ii); delete later
                        ii=ii-1; 
                        break
                    end
                end
            end
            ii=ii+1;
        end
        % delete merged regions
        remove=ireg(merged);
        for ii= remove(end:-1:1)
            cc.PixelIdxList(ii)=[];
        end
    end
  
    
    for nregion=1:length(cc.PixelIdxList)
       
        region=cc.PixelIdxList{nregion};
        

        %% set up east-west equations for weighted inversion
        
        reg=false(1,xi*yi)'; 
        reg(region)=true;
        
        % centered finite differences on a staggered grid
        en= reg & circshift(reg,-yi);
        if strcmp(wrap,'none')  % remove equations for eastern boundary for zonally-nonperiodic domain
            bdy_east=false(1,xi*yi); bdy_east(yi*(xi-1)+1:end)=true;
            en=en & ~bdy_east(:);
        end
               
        % j1 are j-indices for matrix coefficient 1
        sreg=cumsum(reg); % sparse indices of region (points of non-region are indexed with dummy)
        sreg_en=circshift(sreg,-yi); % sparse indices of eastward neighbours
        j1_ew=sreg_en(en); 
       
        j2_ew=sreg(en); % j2 are j-indices for matrix coefficient -1
        

        %% set up north-south equations for weighted inversion
        
        nn= reg & circshift(reg,-1);
        % remove equations for northern boundary
        bdy_north=false(1,xi*yi); bdy_north(yi:yi:yi*xi)=true;
        nn=nn & ~bdy_north(:);

        sreg_nn=circshift(sreg,-1);
        j1_ns=sreg_nn(nn);
        
        j2_ns=sreg(nn);
        
        j2=[j2_ew',j2_ns'];
        
        
        %% make the average of Phi' zero 
        % this should keep the surface from drifting away from the initial condition
        % we might change that to a different condition
        
        j1_condition=[1:sum(reg)];
        j1=[j1_ew',j1_ns',j1_condition];
        
        
        i2=1:(sum(en)+sum(nn)); % i-indices for matrix coeff. -1
        i1=[i2, (sum(en)+sum(nn)+1)*ones(1,sum(reg))];
        
        % build sparse matrices
        A=sparse([i1,i2],[j1,j2],[ones(1,length(i1)),-ones(1,length(i2))]);
        b=sparse( [xx_squeeze(en).*e1t(en); yy_squeeze(nn).*e2t(nn); 0 ]);
        
        disp(['solving for region ',int2str(nregion)]);
        %stef: 'exact' gets to the solution quicker but requires more
        %memory
        switch solver
            case 'iterative'
                [x,dummyflag] = lsqr(A,b,1e-7,50000);
                
            case 'exact'
                x = (A'*A)\(A'*b);
        end
        
        x = full(x)';
        
        % put density changes calculated by the least-squares solver into
        % their appropriate position in the matrix
        
        switch choice
            case 's'
                phiprime(region) = - x;
            case 'epsilon'
                phiprime_e(region) = x;
        end
        
    end

      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dummy_phiprime_e(1,:,:) = (-1).*phiprime_e;
    
    r=1.0;
    delta = 1e-3;

    t2=gsw_rho(sns(:),ctns(:),pns(:));
    t2=t2.*(1+r*dummy_phiprime_e(:)); % rho_s+rho'
    
    inds=1:yi*xi;
    fr=true(1,yi*xi);
    s_tmp=s;
    ct_tmp=ct;
    p_tmp=p;
    pns_tmp = nan(1,yi,xi);
    sns_tmp = nan(1,yi,xi);
    ctns_tmp = nan(1,yi,xi);
    pns=pns(:);

    refine_ints=100;
    
    cnt=0;
    while 1
        cnt=cnt+1;
        inds=inds(fr);
         
        if cnt==1 % in first iteration pns_l is stacked vertically zi times, after that it is stacked refine_ints times
            stack=zi;
        elseif cnt==2
            stack=refine_ints+1;
        end
        if cnt==1 | cnt==2
            ii=bsxfun(@times,1:yi*xi,ones(stack,1)); 
            pns_3d=pns(ii);
            t2_3d_full=t2(ii);
        end
        
        t1=gsw_rho(s_tmp(:,:),ct_tmp(:,:),pns_3d(:,inds));
        t2_3d=t2_3d_full(:,inds);
        tni=t1-t2_3d; % rho-(rho_s+rho')
        
        Itni_p = tni>0;
        Itni_n = tni<0;
        
        zc = any(Itni_p,1) & any(Itni_n,1); % horizontal indices of locations where zero-crossing occurs
        [min_tni, lminr]=min(abs(tni));
        cond1=min_tni>delta; 
        final=min_tni<=delta; 
        fr= zc & cond1; % adjust surface depth here (find root)
        
        lminr=lminr+stack*[0:size(Itni_n,2)-1];
        lminr=lminr(final);
        
        sns_tmp(inds(final))=s_tmp(lminr);
        ctns_tmp(inds(final)) =ct_tmp(lminr);
        pns_tmp(inds(final)) =p_tmp(lminr);


        if all(~fr)
            break      
        end

        k=sum(Itni_n,1);
        k=k+stack*[0:size(Itni_n,2)-1];
        k=k(fr); 
        
        ds_ =  ( s_tmp(k+1) - s_tmp(k))/refine_ints;
        dct_ = (ct_tmp(k+1) - ct_tmp(k))/refine_ints;
        dp_ =  (p_tmp(k+1) - p_tmp(k))/refine_ints;

        ds_ =bsxfun(@times, ds_, [0:refine_ints]');
        dct_ = bsxfun(@times, dct_, [0:refine_ints]');
        dp_ = bsxfun(@times, dp_, [0:refine_ints]');

        s_tmp =  bsxfun(@plus,s_tmp(k),ds_);
        ct_tmp =  bsxfun(@plus,ct_tmp(k),dct_);
        p_tmp =  bsxfun(@plus,p_tmp(k),dp_);

    end
      
    pns=pns_tmp;
    ctns=ctns_tmp;
    sns=sns_tmp;
   
                        
end

    sns_i=sns;
    ctns_i=ctns;
    pns_i=pns;
    
    ithist.sns = sns_hist;
    ithist.ctns = ctns_hist;
    ithist.pns = pns_hist;
    ithist.phiprime_e = phiprime_e_hist;
    ithist.eps_ss = eps_ss;

end



