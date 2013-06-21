function [sns_i,ctns_i,pns_i] = optimize_surface_exact(s,ct,p,g,n2,sns,ctns,pns,e1t,e2t)

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
%   Principal investigator: Trevor McDougall
    %   type 'help analyze_surface' for more information
    %   type 'analyze_surface_license' for license details
    %   type 'analyze_surface_version' for version details
    %

    user_input;
    [zi,yi,xi] = size(s);

    % prepare data
    p_mid = 0.5*(p(2:zi,:,:) + p(1:zi-1,:,:)); % pressure at the grid on which n2 is defined
    cut_off_choice(1,1:yi,1:xi) = mld(s,ct,p); % mixed-layer depth

    % iterations of inversion
    it=0; % start with it=0 and increment it after the initial surface is written to output
    while it<=nit;

        % calculate slope errors/density gradient errors
        [ss,sx,sy,curl_s,ee,ex,ey,curl_e,ver] = slope_error(p,g,n2,sns,ctns,pns,e1t,e2t,'bp'); %#ok

        % diagnose
        if save_iterations;
            if it==0; % dummy value for phiprime
                phiprime=nan;
            end
            diagnose_and_write(it,sns,ctns,pns,ex,ey,phiprime);
        end

        if it==nit; % break out if done
            break
        end

        it=it+1; % start the next iteration
        disp(['Iteration nr.',int2str(it)]);

        % choice between minimizing slope errors or density gradient errors
        switch choice
            case 's'
                xx = sx;
                yy = sy;
            case 'epsilon'
                xx = ex;
                yy = ey;
        end

        % calculate buoyancy frequency on density surface
        n2_ns = var_on_surf(pns,p_mid,n2);

        % disregard data above mixed layer depth
        xx(pns<= cut_off_choice)=nan;
        yy(pns<= cut_off_choice)=nan;
        n2_ns(pns<=cut_off_choice)=nan;
        pns(pns<=cut_off_choice)=nan;

        xx = squeeze(xx);
        yy = squeeze(yy);
        pns = squeeze(pns);
        n2_ns = squeeze(n2_ns);

        % find independent regions -> a least-squares problem is solved for
        % each of these regions
        regions=find_regions(n2_ns);

        % solve for phiprime
        phiprime=solve_lsqr(regions, xx, yy, e1t, e2t);

        % find corrected surface
        dummy_phiprime_e(1,:,:) = (-1).*phiprime;
        if it==1
            r=1+750*dummy_phiprime_e(:);
            r(r<0.9)=0.9; r(r>1.1)=1.1;
        else
            r=1;
        end
        [sns, ctns, pns] = dz_from_phiprime(r,sns, ctns, pns, s, ct, p, phiprime );

    end

    sns_i=sns;
    ctns_i=ctns;
    pns_i=pns;

end


function phiprime=solve_lsqr(regions, xx, yy, e1t, e2t)
    user_input;
    [yi,xi]=size(xx);
    phiprime = nan(yi,xi);

    for iregion=1:length(regions) 
        
        region=regions{iregion};

        % set up east-west equations for weighted inversion
        reg=false(1,xi*yi)'; 
        reg(region)=true;
        
        en= reg & circshift(reg,-yi); % find points where finite differences can be computed
        if ~zonally_periodic;  % remove equations for eastern boundary for zonally-nonperiodic domain
            not_bdy_east=true(1,xi*yi); not_bdy_east(yi*(xi-1)+1:end)=false;
            en=en & not_bdy_east(:);
        end
        sreg=cumsum(reg); % sparse indices of points forming the region (points of non-region are indexed with dummy)
        sreg_en=circshift(sreg,-yi); % sparse indices of eastward neighbours
        
        j1_ew=sreg_en(en);  % j1 are j-indices for matrix coefficient 1
        j2_ew=sreg(en); % j2 are j-indices for matrix coefficient -1
        
        % set up north-south equations for weighted inversion
        nn= reg & circshift(reg,-1);
        not_bdy_north=true(1,xi*yi); not_bdy_north(yi:yi:yi*xi)=false; % remove equations for northern boundary
        nn=nn & not_bdy_north(:);
        sreg_nn=circshift(sreg,-1);
        
        j1_ns=sreg_nn(nn);
        j2_ns=sreg(nn);
           
        % make the average of Phi' zero 
        % this should keep the surface from drifting away from the initial condition
        % we might change that to a different condition
        j1_condition=[1:sum(reg)];
  
        j1=[j1_ew',j1_ns',j1_condition];
        j2=[j2_ew',j2_ns'];
        
        i2=1:(sum(en)+sum(nn)); % i-indices for matrix coeff. -1
        i1=[i2, (sum(en)+sum(nn)+1)*ones(1,sum(reg))]; % i-indices for matrix coeff. 1
        
        % build sparse matrices
        A=sparse([i1,i2],[j1,j2],[ones(1,length(i1)),-ones(1,length(i2))]);
        b=sparse( [xx(en).*e1t(en); yy(nn).*e2t(nn); 0 ]);
        
        disp(['solving for region ',int2str(iregion)]);
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
                phiprime(region) = x;
        end   
    end
end


function [sns_out,ctns_out,pns_out] = dz_from_phiprime(r,sns, ctns, pns, s, ct, p, phiprime_e );

    [zi,yi,xi]=size(s);
    dummy_phiprime_e(1,:,:) = (-1).*phiprime_e;
    
    delta = 1e-9;

    t2=gsw_rho(sns(:),ctns(:),pns(:));
    t2=t2.*(1+r.*dummy_phiprime_e(:)); % rho_s+rho'
    
    inds=1:yi*xi;
    fr=true(1,yi*xi);
    pns_out = nan(1,yi,xi);
    sns_out = nan(1,yi,xi);
    ctns_out = nan(1,yi,xi);
    pns=pns(:);

    refine_ints=100;
    
    cnt=0;
    while 1
        cnt=cnt+1;
         
        if cnt==1 % in first iteration pns_l is stacked vertically zi times, after that it is stacked refine_ints times
            stack=zi;
        elseif cnt==2
            stack=refine_ints+1;
        end
        if cnt==1 | cnt==2
            ii=bsxfun(@times,1:yi*xi,ones(stack,1)); 
            pns_stacked=pns(ii); % stack pressure of current surface vertically
            t2_stacked_full=t2(ii); % stack locally referenced density of current surface vertically
        end
        
        inds=inds(fr); % points where surface has not been corrected
        t1=gsw_rho(s(:,:),ct(:,:),pns_stacked(:,inds)); % 3-d density referenced to pressure of the current surface
        t2_3d=t2_stacked_full(:,inds); %
        tni=t1-t2_3d; % rho-(rho_s+rho'); find corrected surface by finding roots of this term
        
        Itni_p = tni>0;
        Itni_n = tni<0;
        
        zc = any(Itni_p,1) & any(Itni_n,1); % horizontal indices of locations where zero-crossing occurs
        [min_tni, lminr]=min(abs(tni));
        cond1=min_tni>delta; 
        final=min_tni<=delta; % at these horizontal indices root has been found
        fr= zc & cond1; %  at these horizontal locations we have to increase the vertical resolution before finding the root
        
        lminr=lminr+stack*[0:size(Itni_n,2)-1];
        lminr=lminr(final);
        
        sns_out(inds(final))=s(lminr); % adjust surface where root has already been found
        ctns_out(inds(final)) =ct(lminr);
        pns_out(inds(final)) =p(lminr);


        if all(~fr) % break out of loop if all roots have been found
            break      
        end

        k=sum(Itni_n,1); % find indices of flattened 3d-array where vertical resolution must be increased
        k=k+stack*[0:size(Itni_n,2)-1];
        k=k(fr); 
        
        ds_ =  ( s(k+1) - s(k))/refine_ints; % increase resolution in the vertical
        dct_ = (ct(k+1) - ct(k))/refine_ints;
        dp_ =  (p(k+1) - p(k))/refine_ints;

        ds_ =bsxfun(@times, ds_, [0:refine_ints]');
        dct_ = bsxfun(@times, dct_, [0:refine_ints]');
        dp_ = bsxfun(@times, dp_, [0:refine_ints]');

        s =  bsxfun(@plus,s(k),ds_);
        ct =  bsxfun(@plus,ct(k),dct_);
        p =  bsxfun(@plus,p(k),dp_);

    end
end


function regions=find_regions(vv)

    user_input;
    [yi,xi]=size(vv);
    
    % flag wet points (mixed layer excluded)
    wet=~isnan(vv);
    
    cc=bwconncomp(wet,4);

    if zonally_periodic; % in a zonally periodic domain merge regions which are cut apart by grid boundary

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
    
    regions=cc.PixelIdxList;
end


function diagnose_and_write(it,sns,ctns,pns,ex,ey,phiprime_e)
    user_input; % read nit, choice, etc.

    if it==0 % initialize
        [gi,yi,xi]=size(sns);
        slope_square = nan(nit+1,1);

        sns_hist = nan(nit+1,yi,xi); % store variables on initial surface and on nit improvements (=> nit+1)
        ctns_hist = nan(nit+1,yi,xi);
        pns_hist = nan(nit+1,yi,xi);
        eps_rms_hist=nan(nit+1,1);
        phiprime_e_hist = nan(nit,yi,xi);

        
        vars = {'sns_hist','ctns_hist','pns_hist','eps_rms_hist','phiprime_e_hist'};
        save(history_file, vars{:},'-v7.3');
    end
    
    iteration_history = matfile(history_file,'Writable',true);
    iteration_history.sns_hist(it+1,:,:) = sns;
    iteration_history.ctns_hist(it+1,:,:) = ctns;
    iteration_history.pns_hist(it+1,:,:) = pns;   

    if it>0
        iteration_history.phiprime_e_hist(it,:,:) = permute(phiprime_e,[3,1,2]);
    end
    
    if strcmp(choice, 'epsilon')
        square = ex .* ex + ey .* ey; % ex and ey are defined on different grids, but this is not relevant for calculating the root-mean-square
        slope_square(it+1,1) = nansum(square(:));
        no_pts = length(find(~isnan(square(:))));
        iteration_history.eps_rms_hist(it+1,1) = sqrt(slope_square(it+1,1)/no_pts);
    end
    
end



