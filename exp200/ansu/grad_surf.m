function [grad_x,grad_y] = grad_surf(var_surf,e1t,e2t)

%           Find the gradient of a variable in x and y-direction.
%
% Usage:    [grad_x,grad_y] = grad_surf(var_surf,e1t,e2t)
%
%
% Input:    var                   variable 
%           e1t                   zonal scale-factor at var_surf points 
%           e2t                   meridional scale-factor at var_surf 
%                                 points
%
% Output:   grad_x                gradient in zonal direction 
%           grad_y                gradient in meridional direction 
%
% Units:    scale factors         m         
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

user_input; % read switch zonally_periodic
%% check input arguments

if ~(nargin == 3)
    error('grad_surf.m: requires 3 input arguments')
end

%% calculate gradients

[yi,xi] = size(var_surf);

% calculate slope in longitude direction

grad_x = ( circshift(var_surf,[0 -1])-var_surf )./ e1t;
if ~zonally_periodic;
    grad_x(:,xi)=nan;
end

% calculate slope in latitude direction

 grad_y = ( circshift(var_surf,[-1 0])-var_surf )./ e2t;
 grad_y(yi,:) = nan;
 
 
