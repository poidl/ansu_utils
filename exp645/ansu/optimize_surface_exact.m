function [sns_i,ctns_i,pns_i,er_delx,er_dely] = optimize_surface_exact(s,ct,p,sns,ctns,pns)

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

omega_user_input;

%cut_off_choice = mld(s,ct,p); % mixed-layer depth
 
% initialize variables for wetting
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
            er_delx=nan; er_dely=nan; derr=nan; res=nan; b=nan; n2ns=nan; 
        end
        diagnose_and_write(it,sns,ctns,pns,er_delx,er_dely,derr,res,b,n2ns);
    end
    
    %if (it==nit_max) || breakout
    if breakout
        break
    end
    
    % Locations where outcropping occurs may have changed. Add points to
    % surface if necessary.
    apply_wetting;
    
    it=it+1; % start the next iteration
    
    disp(['iter ',int2str(it), '    append ',num2str(nneighbours)]);
        
%     % disregard data above mixed layer depth
%     drhodx(pns<=cut_off_choice)=nan;
%     drhody(pns<=cut_off_choice)=nan;
%     sns(pns<=cut_off_choice)=nan;
%     ctns(pns<=cut_off_choice)=nan;
%     pns(pns<=cut_off_choice)=nan;

    if no_land_mask
        load([datapath,'latlon.mat'])
        [ocean, n] = gamma_ocean_and_n(s,ct,p,lon,lat);
        save([datapath,'no_land_mask.mat'], 'ocean', 'n')
    end
 
    if strcmp(error_measure,'drho_local')
        [er_delx,er_dely]=delta_tilde_rho(sns,ctns,pns);
    elseif strcmp(error_measure,'slope_difference')
        [er_delx,er_dely]=delta_slope(sns,ctns,pns,s,ct,p);
    end
%     %if strcmp(error_measure,'drho_local')
%         [er_delx,er_dely]=delta_tilde_rho(sns,ctns,pns);
%     %elseif strcmp(error_measure,'slope_difference')
%         [er_delx2,er_dely2]=delta_slope(sns,ctns,pns,s,ct,p);
%     %end
    
    
    if strcmp(error_measure,'drho_local')
        
        if use_b
            error('no area weighting')
            [er_delx,er_dely,regions,b]=use_bstar(er_delx,er_dely,pns,s,ct,p);
            [derr,res]=solve_lsqr(regions, er_delx, er_dely);
            derr=derr./b;
            
        else
            if strcmp(weights,'drho_local_n2')
                if it==5
                    [n2ns,n2nsx,n2nsy]=get_n2ns(pns,s,ct,p);
                end
                if it>=5
                    erx=erx./n2nsx;
                    ery=ery./n2nsy;
                end
                
            end
        
            % We want to minimize an area integral. Take account of grid
            % point spacing.
            load([datapath,'dxdy.mat']) 
            er_gradx=er_delx./dx;
            er_grady=er_dely./dy;
            [erx,ery]=times_sqrtdA_on_delta(er_gradx,er_grady,dx,dy);             


            regions=find_regions(pns);
            [derr,res]=solve_lsqr(regions, erx, ery);
            derr=0.8*derr;
            
            if strcmp(weights,'drho_local_n2') & it>=5
                derr=0.2*derr.*n2ns;
            end    
        end % use_b else
        [sns, ctns, pns] = depth_ntp_simple(sns(:)', ctns(:)', pns(:)', s(:,:), ct(:,:), p(:,:), derr(:)' );
        [zi,yi,xi]=size(s);
        sns=reshape(sns,[yi xi]);
        ctns=reshape(ctns,[yi xi]);
        pns=reshape(pns,[yi xi]);
        
    elseif strcmp(error_measure,'slope_difference') 
        
        % We want to minimize an area integral. Take account of grid
        % point spacing.
        load([datapath,'dxdy.mat']) 
        er_gradx=er_delx./dx;
        er_grady=er_dely./dy;
        [erx,ery]=times_sqrtdA_on_delta(er_gradx,er_grady,dx,dy);        
        
        [regions]=remove_points(erx,ery,pns); % not possible to use find_regions() here, since the differencing may have introduced nans
        [derr,res]=solve_lsqr(regions, erx, ery);
       
        pns=pns+0.1*derr;
        sns=var_on_surf_stef(s,p,pns);
        ctns=var_on_surf_stef(ct,p,pns);
        
    end
       
    % remove any regions which may have been detached during optimization
    %keyboard
    if clamp_on_point
        load([datapath,'ibb.mat'])
        [sns,ctns,pns] = get_connected(sns,ctns,pns,ibb);
    end

end

sns_i=sns;
ctns_i=ctns;
pns_i=pns;

end
