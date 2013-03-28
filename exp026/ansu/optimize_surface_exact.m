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

switch choice
    case 's'
        keyboard
    case 'epsilon'
        [ex_o] = cut_off(ex,pns,cut_off_choice);
        [ey_o] = cut_off(ey,pns,cut_off_choice);
        square_o = ex_o .* ex_o + ey_o .* ey_o;
        slope_square_o = nansum(square_o(:));
        no_pts = length(find(~isnan(square_o(:))));
        % eps_ss(1,1) = sqrt(slope_square_o/no_pts);
end

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
depth_change_i_hist = nan(nit,yi,xi);
r_hist = nan(nit);
eps_ss=nan(nit,1);
ppm = nan(nit,4);

%% iterations of inversion
it = 0;
phiprime = 0;


for it = 1:nit
    
    %it = it + 1;
    
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
    
    
    b = zeros(3*xi*yi,1,1);
    s1 = zeros(3*xi*yi,1);
    s2 = zeros(3*xi*yi,1);
    s3 = zeros(3*xi*yi,1);
    
    % find grid points which are not continents
    
    inds = find(~isnan(n2_ns_squeeze(:)));
    
    % build matrix where the ocean gridpoints have their indices and
    % continet gridpoints nans
    
    ng = nan(yi,xi);
    ng(inds) = inds;
    
    % select regions
    
    % build neighbour matrix -> build a matrix which is like a look-up
    % table to see which gridpoints communicate with each other
    
    neighbour = nan(5,length(inds));
    
    switch wrap
        
        case 'none'
            
            for i = 1:length(inds)
                
                [jj,ii] = ind2sub([yi,xi],inds(i));
                
                if (jj == 1) && (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = nan;
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (jj == yi) && (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = nan;
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == 1) && (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = nan;
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (jj == yi) && (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = nan;
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = nan;
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == yi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = nan;
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                else
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                end
            end
            
        case 'long'
            
            for i = 1:length(inds)
                
                [jj,ii] = ind2sub([yi,xi],inds(i));
                
                if (jj == 1) && (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,xi);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (jj == yi) && (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == 1) && (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (jj == yi) && (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,xi);
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = nan;
                elseif (ii == 1)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,xi);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (jj == yi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = nan;
                    neighbour(5,i) = ng(jj-1,ii);
                elseif (ii == xi)
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                else
                    neighbour(1,i) = ng(jj,ii);
                    neighbour(2,i) = ng(jj,ii+1);
                    neighbour(3,i) = ng(jj,ii-1);
                    neighbour(4,i) = ng(jj+1,ii);
                    neighbour(5,i) = ng(jj-1,ii);
                end
            end
    end
    
    nregion = 0;
    region_matrix = nan(yi,xi);
    
    % find independent regions -> a least-squares problem is solved for
    % each of these regions
    
    while (length(find(~isnan(region_matrix))) ~= length(neighbour))
        
        pos = 0;
        %region = [];
        region = nan(1,2*length(inds));
        nregion = nregion + 1;
        
        % find starting point of region
        
        for i = 1:length(inds)
            
            if ~isnan(neighbour(1,i)) && isnan(region_matrix(neighbour(1,i)))
                pos = pos+1;
                region(pos) = neighbour(1,i);
                region_matrix(neighbour(1,i)) = nregion;
                neighbour(1,i) = nan;
                
                if ~isnan(neighbour(2,i))
                    if isnan(region_matrix(neighbour(2,i)))
                        pos = pos+1;
                        region(pos) = neighbour(2,i);
                        region_matrix(neighbour(2,i)) = nregion;
                        neighbour(2,i) = nan;
                    end
                end
                
                if ~isnan(neighbour(3,i))
                    if isnan(region_matrix(neighbour(3,i)))
                        pos = pos+1;
                        region(pos) = neighbour(3,i);
                        region_matrix(neighbour(3,i)) = nregion;
                        neighbour(3,i) = nan;
                    end
                end
                
                if ~isnan(neighbour(4,i))
                    if isnan(region_matrix(neighbour(4,i)))
                        pos = pos+1;
                        region(pos) = neighbour(4,i);
                        region_matrix(neighbour(4,i)) = nregion;
                        neighbour(4,i) = nan;
                    end
                end
                
                if ~isnan(neighbour(5,i))
                    if isnan(region_matrix(neighbour(5,i)))
                        pos = pos+1;
                        region(pos) = neighbour(5,i);
                        region_matrix(neighbour(5,i)) = nregion;
                        neighbour(5,i) = nan;
                    end
                end
                
                break
                
            end
        end
        
        region_old = [];
        
        %while(length(region) ~= length(region_old))
        while(length(region(1:pos)) ~= length(region_old))
            
            region_old = region(1:pos);
            
            % find all other points of region
            
            for i = 1:length(inds)
                
                if ~isnan(neighbour(1,i))
                    
                    if (region_matrix(neighbour(1,i)) == nregion) %#ok
                        
                        if ~isnan(neighbour(1,i))
                            if isnan(region_matrix(neighbour(1,i)))
                                pos = pos+1;
                                region(pos) = neighbour(1,i);
                                region_matrix(neighbour(1,i)) = nregion;
                                neighbour(1,i) = nan;
                            end
                        end
                        
                        if ~isnan(neighbour(2,i))
                            if isnan(region_matrix(neighbour(2,i)))
                                pos = pos+1;
                                region(pos) = neighbour(2,i);
                                region_matrix(neighbour(2,i)) = nregion;
                                neighbour(2,i) = nan;
                            end
                        end
                        
                        if ~isnan(neighbour(3,i))
                            if isnan(region_matrix(neighbour(3,i)))
                                pos = pos+1;
                                region(pos) = neighbour(3,i);
                                region_matrix(neighbour(3,i)) = nregion;
                                neighbour(3,i) = nan;
                            end
                        end
                        
                        if ~isnan(neighbour(4,i))
                            if isnan(region_matrix(neighbour(4,i)))
                                pos = pos+1;
                                region(pos) = neighbour(4,i);
                                region_matrix(neighbour(4,i)) = nregion;
                                neighbour(4,i) = nan;
                            end
                        end
                        
                        if ~isnan(neighbour(5,i))
                            if isnan(region_matrix(neighbour(5,i)))
                                pos = pos+1;
                                region(pos) = neighbour(5,i);
                                region_matrix(neighbour(5,i)) = nregion;
                                neighbour(5,i) = nan;
                            end
                        end
                    end
                end
            end
        end
        region(pos+1:length(region)) = [];
        % only use points of the region improved in this loop
        
        xx_inds = xx_squeeze(region);
        yy_inds = yy_squeeze(region);
        e1t_inds = e1t(region);
        e2t_inds = e2t(region);
        
        % build a matrix where all the points of the region are labelled
        % with 1,2,....length(region), other regions and continents are
        % filled with nans
        
        ng_region = nan(yi,xi);
        ng_region(region) = 1:length(region);
        
        %% set up east-west equations for weighted inversion
        
        neq = 0; % set equation number to 0
        ieq = 0;
        
        for i = 1:length(region)
            
            [jj,ii] = ind2sub([yi,xi],region(i));
            
            switch wrap
                
                case {'none'}
                    
                    if (ii+1 <= xi) && (~isempty(find(region == ng(jj,ii)))) && (~isempty(find(region == ng(jj,ii+1)))) && (~isnan(xx_inds(ng_region(jj,ii)))) %#ok
                        
                        neq = neq + 1;
                        ieq = ieq + 1;
                        s1(ieq) = neq;
                        s2(ieq) = ng_region(jj,ii);
                        s3(ieq) = -1;
                        ieq = ieq+1;
                        s1(ieq) = neq;
                        s2(ieq) = ng_region(jj,ii+1);
                        s3(ieq) = 1;
                        b(neq,1) = xx_inds(ng_region(jj,ii)) * e1t_inds(ng_region(jj,ii));
                        
                    end
                    
                case {'long'}
                    
                    if (ii+1 <= xi)
                        
                        if (~isempty(find(region == ng(jj,ii)))) && (~isempty(find(region == ng(jj,ii+1)))) && (~isnan(xx_inds(ng_region(jj,ii)))) %#ok
                            
                            neq = neq + 1;
                            ieq = ieq + 1;
                            s1(ieq) = neq;
                            s2(ieq) = ng_region(jj,ii);
                            s3(ieq) = -1;
                            ieq = ieq+1;
                            s1(ieq) = neq;
                            s2(ieq) = ng_region(jj,ii+1);
                            s3(ieq) = 1;
                            b(neq,1) = xx_inds(ng_region(jj,ii)) * e1t_inds(ng_region(jj,ii));
                            
                        end
                        
                    elseif (ii+1 > xi)
                        
                        if (~isempty(find(region == ng(jj,ii))) && ~isempty(find(region == ng(jj,1)))) && (~isnan(xx_inds(ng_region(jj,ii))))%#ok
                            
                            neq = neq + 1;
                            ieq = ieq + 1;
                            s1(ieq) = neq;
                            s2(ieq) = ng_region(jj,ii);
                            s3(ieq) = -1;
                            ieq = ieq+1;
                            s1(ieq) = neq;
                            s2(ieq) = ng_region(jj,1);
                            s3(ieq) = 1;
                            b(neq,1) = xx_inds(ng_region(jj,ii)) * e1t_inds(ng_region(jj,ii));
                            
                        end
                        
                    end
            end
        end
        
        %% set up north-south equations for weighted inversion
        
        for i = 1:length(region)
            
            [jj,ii] = ind2sub([yi,xi],region(i));
            
            if (jj+1 <= yi) && (~isempty(find(region == ng(jj,ii)))) && ~isempty(find(region == ng(jj+1,ii))) && (~isnan(yy_inds(ng_region(jj,ii))))%#ok
                
                neq = neq + 1;
                ieq = ieq + 1;
                s1(ieq) = neq;
                s2(ieq) = ng_region(jj,ii);
                s3(ieq) = -1;
                ieq = ieq + 1;
                s1(ieq) = neq;
                s2(ieq) = ng_region(jj+1,ii);
                s3(ieq) = 1;
                b(neq,1) = yy_inds(ng_region(jj,ii)) * e2t_inds(ng_region(jj,ii));
            end
        end
        
        % make the average of all density changes zero -> this should keep
        % the surface from drifting away from the initial condition
        
        neq = neq + 1;
        %stef: this can stop the surface of drifting away. We might change
        % that to a different condition.
        for i = 1:length(region)
            [jj,ii] = ind2sub([yi,xi],region(i));
            ieq = ieq + 1;
            s1(ieq) = neq;
            s2(ieq) = ng_region(jj,ii);
            s3(ieq) = 1;
        end
        
        b(neq,1) = 0;
        
        % cut A and b to appropriate size
        
        s1 = s1(1:ieq);
        s2 = s2(1:ieq);
        s3 = s3(1:ieq);
        b = b(1:neq,1);
        
        % make matrix sparse and invert
        
        A = sparse(s1,s2,s3);
        b = sparse(b);
        
        % disp(['solving for region ',int2str(nregion)]);
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
    r=0.5;
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



