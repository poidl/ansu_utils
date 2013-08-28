function [ex,ey] = epsilon(p,sns,ctns,pns,e1t,e2t)

%           Calculate density gradient errors
%
% Usage:    [ex,ey] = epsilon(p,sns,ctns,pns,e1t,e2t)
%
%           Calculate slope errors, density gradient errors, their curl and
%           the fictitious diapycnal diffusivity
%
% Input:    p           pressure
%           sns         salinity on neutral surface
%           ctns        conservative temperature on neutral surface
%           pns         pressure on neutral surface
%           e1t         zonal scale-factor at tracer points
%           e2t         meridional scale-factor at tracer points
%
% Output:  
%           ex         x-component of epsilon
%           ey         y-component of epsilon
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

%% check input arguments

if ~(nargin == 6)
    error('epsilon.m: requires 6 input arguments')
end

%% initialize

[zi,yi,xi] = size(p);

%% calculate slope errors and epsilon

[gradx_ct,grady_ct] = grad_surf(ctns,e1t,e2t);
[gradx_s,grady_s] = grad_surf(sns,e1t,e2t);

[yy,xx]=size(sns);
[dummy,alpha,beta]=gsw_rho_alpha_beta(sns(:,:),ctns(:,:),pns(:,:));
rho=gsw_rho(sns(:,:),ctns(:,:),pns(:,:));
alpha=rho.*alpha;
beta=rho.*beta;


alpha=reshape(alpha,[yy,xx]);
beta=reshape(beta,[yy,xx]);

alpha_x = 0.5*( alpha + circshift(alpha,[0 -1]) );
alpha_y = 0.5*( alpha + circshift(alpha,[-1 0]) );
beta_x = 0.5 * ( beta + circshift(beta,[0 -1]) );
beta_y = 0.5 * ( beta + circshift(beta,[-1 0]) );

alpha_y(yi,:) = nan;
beta_y(yi,:) = nan;

if ~zonally_periodic;
    alpha_x(:,xi) = nan;
    beta_x(:,xi) = nan;
end

% calculate density gradient errors (epsilon)

ex = ((beta_x .* gradx_s) - (alpha_x .* gradx_ct));
ey = ((beta_y .* grady_s) - (alpha_y .* grady_ct));

