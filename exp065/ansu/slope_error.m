function [ss,sx,sy,curl_s,ee,ex,ey,curl_e,fdd] = slope_error(p,g,n2,sns,ctns,pns,e1t,e2t,keyword)

%           Calculate slope errors
%
% Usage:    [ss,sx,sy,curl_s,ee,ex,ey,curl_e,ver] = slope_error(p,g,n2,sns,ctns,pns,e1t,e2t,keyword)
%
%           Calculate slope errors, density gradient errors, their curl and
%           the fictitious diapycnal diffusivity
%
% Input:    p           pressure
%           g           gravitational accelleration
%           sns         salinity on neutral surface
%           ctns        conservative temperature on neutral surface
%           pns         pressure on neutral surface
%           n2_ns       buoyancy frequency
%           e1t         zonal scale-factor at tracer points
%           e2t         meridional scale-factor at tracer points
%           lats        latitude
%           keyword     'op' for values on points
%                       'bp' for values in between points
%
% Output:   ss         slope error
%           sx         x-component of slope error
%           sy         y-component of slope error
%           curl_s     curl of slope error
%           ee         epsilon
%           ex         x-component of epsilon
%           ey         y-component of epsilon
%           curl_e     curl of epsilon
%           fdd        fictitious diapycnal diffusivity
%
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%           buoyancy frequency          s^-1
%           gravitational acceleration  m/s^2
%           scale factors               m
%           fictitious diapycnal diff.  m^2/s
%           
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   Partially modified by P. Barker (2010-13)
%   Partially modified by S. Riha (2013)
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

user_input;

%% check input arguments

if ~(nargin == 9)
    error('slope_error.m: requires 9 input arguments')
end

%% initialize

[gi,dummy,dummy] = size(sns); %#ok
[zi,yi,xi] = size(p);

%% make g correct size

if (length(size(g)) == 2)
    g=permute(g,[3 1 2]);
end

%% calculate slope errors and epsilon

if strcmp(keyword,'op')

    error('not implemented')
    
elseif strcmp(keyword,'bp')

        % ss, curl_s, ee, curl_e and fictitious diapycnal diffusivity are not calculated due to sx, sy, ex, ey being defined on different grids;

        [gradx_ct,grady_ct] = grad_surf(ctns,e1t,e2t,'bp');
        [gradx_s,grady_s] = grad_surf(sns,e1t,e2t,'bp');
        
        [zz,yy,xx]=size(sns);
        [dummy,alpha,beta]=gsw_rho_alpha_beta(sns(:,:),ctns(:,:),pns(:,:));
        alpha=reshape(alpha,[zz,yy,xx]);
        beta=reshape(beta,[zz,yy,xx]);

        alpha_x = 0.5*( alpha + circshift(alpha,[0 0 -1]) );
        alpha_y = 0.5*( alpha + circshift(alpha,[0 -1 0]) );
        beta_x = 0.5 * ( beta + circshift(beta,[0 0 -1]) );
        beta_y = 0.5 * ( beta + circshift(beta,[0 -1 0]) );
       
        p_mid = (p(2:zi,:,:) + p(1:zi-1,:,:)) ./ 2;
        n2_ns = var_on_surf(pns,p_mid,n2);
        n2_ns(n2_ns==0)=nan;
        fac_tmp = g ./ n2_ns;

        fac_x = 0.5 * ( fac_tmp + circshift(fac_tmp,[0 0 -1]) );
        fac_y = 0.5 * ( fac_tmp + circshift(fac_tmp,[0 -1 0]) );
        
        alpha_y(:,yi,:) = nan;
        beta_y(:,yi,:) = nan;
        fac_y(:,yi,:) = nan;

        if ~zonally_periodic;

            alpha_x(:,:,xi) = nan;
            beta_x(:,:,xi) = nan;
            fac_x(:,:,xi) = nan;
            
        end

        % calculate density gradient errors (epsilon)

        ex = ((beta_x .* gradx_s) - (alpha_x .* gradx_ct));
        ey = ((beta_y .* grady_s) - (alpha_y .* grady_ct));
        ee = nan; % ex and ey are at different locations

        % calculate slope errors

        sx = fac_x .* ex;
        sy = fac_y .* ey;
        ss = nan; % sx and sy are at different locations

        fdd = nan;

else
        error('slope_error.m: keyword must be "op" or "bp"')
end

%% calculate curl of slope errors and epsilon

if strcmp(keyword,'op')

        error('not implemented')

else

        curl_e = [];
        curl_s = [];
end