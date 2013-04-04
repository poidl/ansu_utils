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
% Calls:    cut_off.m, mld.m, slope_error.m, var_on_surf.m
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

%% get size of 3-dim hydrography
%global settings lats longs
[zi,yi,xi] = size(s);

%% calculate buoyancy frequency on density surface

p_mid = 0.5*(p(2:zi,:,:) + p(1:zi-1,:,:));
n2_ns = var_on_surf(pns,p_mid,n2);

%% calculate slope errors/density gradient errors on initial density surface

[ss,sx,sy,curl_s,ee,ex,ey,curl_e,ver] = slope_error(p,g,n2,sns,ctns,pns,e1t,e2t,'bp',wrap); %#ok

%% choice between minimizing slope errors or density gradient errors

switch choice
    case 's'
        xx = sx;
        yy = sy;
    case 'epsilon'
        xx = ex;
        yy = ey;
end

%% disregard data above mixed layer depth

cut_off_choice(1,1:yi,1:xi) = mld(s,ct,p);

[pns] = cut_off(pns,pns,cut_off_choice);
[n2_ns] = cut_off(n2_ns,pns,cut_off_choice);


%% prepare data

slope_square = nan(nit,1);
pns_squeeze = squeeze(pns);
n2_ns_squeeze = squeeze(n2_ns);
xx_squeeze = squeeze(xx);
yy_squeeze = squeeze(yy);

sns_i_hist = nan(nit,yi,xi);
ctns_i_hist = nan(nit,yi,xi);
pns_i_hist = nan(nit,yi,xi);
depth_change_e_i_hist = nan(nit,yi,xi);
eps_ss=nan(nit,1);
ppm = nan(nit,4);

%% iterations of inversion

for it = 1:nit
    
    disp(['Iteration nr.',int2str(it)]); % print number of iteration
    
    if (it > 1)
        
        % disregard data above mixed layer depth
        
        [pns_i] = cut_off(pns_i,pns_i,cut_off_choice);
        [n2_ns_i] = cut_off(n2_ns_i,pns_i,cut_off_choice);
        
        % prepare data for next iteration
        
        pns_squeeze = squeeze(pns_i);
        n2_ns_squeeze = squeeze(n2_ns_i);
        
        switch choice
            case 'epsilon'
                [ex_i] = cut_off(ex_i,pns_i,cut_off_choice);
                [ey_i] = cut_off(ey_i,pns_i,cut_off_choice);
                xx_squeeze = squeeze(ex_i);
                yy_squeeze = squeeze(ey_i);
            case 's'
                [sx_i] = cut_off(sx_i,pns_i,cut_off_choice);
                [sy_i] = cut_off(sy_i,pns_i,cut_off_choice);
                xx_squeeze = squeeze(sx_i);
                yy_squeeze = squeeze(sy_i);
        end
    end
    
    [yi,xi] = size(pns_squeeze); % find dimensions in lats and longs
    
    % preallocate memory
    
    depth_change = nan(yi,xi);
    depth_change_e = nan(yi,xi);
    
    % flag wet points
    
    wet=~isnan(n2_ns_squeeze);
    
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
                    % check if western boarder of region iw(ii) intersects
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
        
        % for forward finite differences only keep equations for points 
        % wich are in the region and have an eastward neighbour in the region 
        en= reg & circshift(reg,-yi);
        if strcmp(wrap,'none')  % remove equations for eastern boundary for zonally-nonperiodic domain
            bdy_east=false(1,xi*yi); bdy_east(yi*(xi-1)+1:end)=true;
            en=en & ~bdy_east(:);
        end
               
        % j1 are j-indices for matrix coefficient 1
        sreg=cumsum(reg); % sparse indices of region (points of non-region are indexed with dummy)
        sreg_en=circshift(sreg,-yi); % sparse indices of eastward neighbours
        j1_ew=sreg_en(en); % keep only those with eastward neighbour in the region
       
        j2_ew=sreg(en); % j2 are j-indices for matrix coefficient -1
        

        %% set up north-south equations for weighted inversion
        
        nn= reg & circshift(reg,-1); % true if northward neighbour is in the region
        % remove equations for northern boundary
        bdy_north=false(1,xi*yi); bdy_north(yi:yi:yi*xi)=true;
        nn=nn & ~bdy_north(:);

        sreg_nn=circshift(sreg,-1);
        j1_ns=sreg_nn(nn);
        
        j2_ns=sreg(nn);
        
        j2=[j2_ew',j2_ns'];
        
        
        %% make the average of all density changes zero 
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
                depth_change(region) = - x;
            case 'epsilon'
                depth_change_e(region) = x;
        end
        
        eval(['region_',int2str(nregion),' = region;']);
        
    end
    
%     if (it > 1)
%         pns_old = pns_i;
%     else
%         pns_old = pns;
%     end
    
    if (it > 1)
        pns_l = pns_i;
        sns_l = sns_i;
        ctns_l = ctns_i;
    else
        pns_l = pns;
        sns_l = sns;
        ctns_l = ctns;
    end

    % calculate new depth of approximate neutral surface
    
    %keyboard
    [zi,yi,xi] = size(s);
    %[ms1 ms2 ms3] = size(s);
    clear dummy_depth_change_e
    dummy_depth_change_e(1,:,:) = (-1).*depth_change_e;
    clear tni
    tmp1=nan*ones(size(s));
    for kk=1:size(s,1)
        tmp1(kk,:,:)=gsw_rho(squeeze(s(kk,:,:)),squeeze(ct(kk,:,:)),squeeze(pns_l));
    end
    tmp2=gsw_rho(squeeze(sns_l),squeeze(ctns_l),squeeze(pns_l));
    tmp3=repmat(permute(tmp2,[3,1,2]),[zi,1,1]);
    r=1.0;
    tni= tmp1- tmp3.*(1+r*repmat(dummy_depth_change_e, [zi,1,1])); % rho-(rho_s+rho')
    
    delta = 1e-6;
    
    pns_i = nan(1,yi,xi);
    ctns_i = nan(1,yi,xi);
    sns_i = nan(1,yi,xi);
    %stef: this finds the dz from qhi'
    for ii_tni = 1:yi
        for jj_tni = 1:xi
            [Itni] = find(~isnan(tni(:,ii_tni,jj_tni)));
            if ~isempty(Itni)
                tni_temp = tni(Itni,ii_tni,jj_tni);
                s_temp = s(Itni,ii_tni,jj_tni);
                ct_temp =ct(Itni,ii_tni,jj_tni);
                p_temp = p(Itni,ii_tni,jj_tni);
                tnif = 0;
                tni_count = 0;
                while tnif == 0
                    tni_count = tni_count + 1;
                    if min(abs(tni_temp)) > delta % min(rho-rho_s-rho') 
                        Itni_p = find(tni_temp > 0);
                        Itni_n = find(tni_temp < 0);
                        if ~isempty(Itni_p) & ~isempty(Itni_n)
                            [dummy tni_ui] = max(tni_temp(Itni_n)); %stef: maximum of negative tni_temp ui 'upper index'
                            [dummy tni_li] = min(tni_temp(Itni_p));   %stef: minimum of positive tni_tmp li 'lower index'
                            if tni_count > 100
                                tni_li = tni_ui + 10;
                            end
                            %stef: devide by 100 instead of iteration, which is slower
                            ii1=Itni_p(tni_li);
                            ii2=Itni_n(tni_ui);
                            ds_ =  ( s_temp(ii1) - s_temp(ii2))/100;
                            dct_ = (ct_temp(ii1) - ct_temp(ii2))/100;
                            dp_ =  (p_temp(ii1) - p_temp(ii2))/100;
                            s_dummy =  [s_temp(ii2):    ds_   :s_temp(ii1)];
                            ct_dummy = [ct_temp(ii2):  dct_  :ct_temp(ii1)];
                            p_dummy =  [p_temp(ii2):   dp_  :p_temp(ii1)];
                            
                            if isempty(s_dummy) | isempty(ct_dummy) | isempty(p_dummy)
                                expand_tni = 0;
                                try
                                    while expand_tni == 0
                                        tni_li = tni_li + 1;
                                        ii1=Itni_p(tni_li);
                                        ds_ =  ( s_temp(ii1) - s_temp(ii2))/100;
                                        dct_ = (ct_temp(ii1) - ct_temp(ii2))/100;
                                        dp_ =  (p_temp(ii1) - p_temp(ii2))/100;
                                        s_dummy = [s_temp(ii2): ds_ :s_temp(ii1)];
                                        ct_dummy = [ct_temp(ii2): dct_ :ct_temp(ii1)];
                                        p_dummy = [p_temp(ii2): dp_ :p_temp(ii1)];
                                        if ~isempty(s_dummy) & ~isempty(ct_dummy) & ~isempty(p_dummy)
                                            expand_tni = 1;
                                        end
                                    end
                                catch
                                    pns_i(1,ii_tni,jj_tni) = NaN;
                                    tnif = 1;
                                    expand_tni = 1;
                                end
                            end
                            if tnif ~= 1
                                ms1 = length(s_dummy);
                                
                                tmp1=gsw_rho(s_dummy(:),ct_dummy(:),repmat(pns_l(1,ii_tni,jj_tni),[ms1,1,1]));
                                tmp2=repmat(gsw_rho(sns_l(1,ii_tni,jj_tni),ctns_l(1,ii_tni,jj_tni),pns_l(1,ii_tni,jj_tni)),[ms1,1,1]);
                                tni_temp=tmp1-tmp2.*(1+r*repmat(dummy_depth_change_e(1,ii_tni,jj_tni),[ms1,1,1]));
                                                                
                                s_temp = s_dummy;
                                ct_temp = ct_dummy;
                                p_temp = p_dummy;
                                
                            end
                        else
                            pns_i(1,ii_tni,jj_tni) = NaN;
                            tnif = 1;
                        end
                    else
                        %stef: save fields from every iteration
                        [dummy Iminr] = min(abs(tni_temp));
                        pns_i(1,ii_tni,jj_tni) = p_temp(Iminr);
                        ctns_i(1,ii_tni,jj_tni) = ct_temp(Iminr);
                        sns_i(1,ii_tni,jj_tni) = s_temp(Iminr);
                        
                        tnif = 1;
                    end
                end
            end
        end
    end
      
    % calculate slope errors and buoyancy frequency on new
    % approximate neutral surface
    
    n2_ns_i = var_on_surf(pns_i,p_mid,n2);
    
    [ss_i,sx_i,sy_i,curl_s_i,ee_i,ex_i,ey_i,curl_e_i] = ...
        slope_error(p,g,n2,sns_i,ctns_i,pns_i,e1t,e2t,'bp',wrap); %#ok
    
    
    % calculate epsilon^2 to estimate quality of approximate neutral
    % surface
    
    switch choice
        case 'epsilon'
            square = ex_i .* ex_i + ey_i .* ey_i;
            slope_square(it,1) = nansum(square(:));
            no_pts = length(find(~isnan(square(:))));
            eps_ss(it,1) = sqrt(slope_square(it,1)/no_pts);
            
        case 's'
            square = sx_i .* sx_i + sy_i .* sy_i;
            slope_square(it,1) = nansum(square(:));
    end
   
    % keyboard
     
    sns_i_hist(it,:,:) = squeeze(sns_i);
    ctns_i_hist(it,:,:) = squeeze(ctns_i);
    pns_i_hist(it,:,:) = squeeze(pns_i(1,:,:));
    depth_change_e_i_hist(it,:,:) = squeeze(depth_change_e);
    jnk1 = (depth_change_e(:));
    ppm(it,1) = mean(jnk1(find(~isnan(jnk1)))) - 2*std(jnk1(find(~isnan(jnk1))));
    ppm(it,2) = mean(jnk1(find(~isnan(jnk1)))) + 2*std(jnk1(find(~isnan(jnk1))));
    ppm(it,3) = std(abs(jnk1(find(~isnan(jnk1)))));
    ppm(it,4) = nansum(abs(jnk1(find(~isnan(jnk1)))))/length(find(~isnan(jnk1)));
    phi_prime_rms = ppm(it,2) - ppm(it,1);
                        
end

    ithist.sns_i = sns_i_hist;
    ithist.ctns_i = ctns_i_hist;
    ithist.pns_i = pns_i_hist;
    ithist.depth_change_e_i = depth_change_e_i_hist;
    ithist.eps_ss = eps_ss;

end



