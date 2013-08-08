function [mld_dpth] = mld(s,ct,p)

%           Calculate mixed-layer depth
%
% Usage:    [mld_dpth] = mld(s,ct,p)
%
%           Calculate the mixed-layer depth as in 'Mixed layer depth over the 
%           global ocean: An examination of profile data and a profile-based 
%           climatology, de Boyer Montégut et al., JGR, 2004'. The criterion
%           selected is a threshold value of density from a surface value 
%           (Δrho = 0.03 kg m−3). 
%
% Input:    s         salinity 
%           ct        conservative temperature 
%           p         pressure 
%
% Output:   mld       mixed-layer depth 
%  
% Calls:    gsw_rho_CT_exact.m
%
% Units:    salinity                 psu (IPSS-78)
%           conservative temperature degrees C (IPS-90)
%           pressure                 dbar
%           mld                      dbar
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


%% check input arguments

if ~(nargin == 3)
  error('mld.m: requires 3 input arguments')
end

[zi,yi,xi] = size(s);

%% calculate gamma^rf

dens = gsw_rho(s(:,:),ct(:,:),zeros(zi,yi*xi));

%% calculate mixed layer depth

min_dens=dens(1,:);
mld_dpth=nan*ones(1,yi*xi);
thresh=bsxfun(@times,min_dens+0.3,ones(zi,1));
pos=(thresh-dens)>0;

cs=cumsum(flipud(pos),1)+1; 
cs(cs~=1)=0;
ip=zi-sum(cs,1);
% replacing the above three lines by
% ip=sum(pos,1);
% will give the wrong answer for gk_ak_gamma.mat
% I think this is because 'pos' can be false on top and turn true deeper. better to start counting falses from bottom.

ii1=ip+zi*[0:(yi*xi)-1];
mld_dpth(ip>0)=p(ii1(ip>0));
mld_dpth=reshape(mld_dpth,[yi,xi]);



