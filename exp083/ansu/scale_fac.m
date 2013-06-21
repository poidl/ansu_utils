function [e2t,e1t] = scale_fac(lats,longs)

%           Find distances between gridpoints of given latitude/longitude 
%
% Usage:    [e2t,e1t] = scale_fac(lats,longs)
%
%           Find distances between gridpoints of given latitude/longitude 
%           in a 2-dim dataset 
%
% Input:    lats        latitude
%           longs       longitude
%
% Output:   e2t         distance between gridpoints in N-S direction             
%           e1t         distance between gridpoints in E-W direction 
%
% Units:    e2t         m
%           e1t         m
%
% Calls:    gsw_distance.m
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

if ~(nargin == 2)
    error('scale_fac.m: requires 2 input arguments')
end

%% initialize and preallocate memory

[yi,xi] = size(lats);
e1t = nan(yi,xi);
e2t = nan(yi,xi);

%% calculate distances

if ~zonally_periodic;
        for j = 1:yi
            for i = 1:xi-1
                e1t(j,i) = gsw_distance([longs(j,i) longs(j,i+1)],[lats(j,i) lats(j,i)]);
            end
        end
        for j = 1:yi-1
            for i = 1:xi
                e2t(j,i) = gsw_distance([longs(j,i) longs(j,i)],[lats(j,i) lats(j+1,i)]);
            end
        end
else
        for j = 1:yi
            for i = 1:xi-1
                e1t(j,i) = gsw_distance([longs(j,i) longs(j,i+1)],[lats(j,i) lats(j,i)]);
            end
        end
        for j = 1:yi
            e1t(j,xi) = gsw_distance([longs(j,xi) longs(j,1)],[lats(j,i) lats(j,i)]);
        end
        for j = 1:yi-1
            for i = 1:xi
                e2t(j,i) = gsw_distance([longs(j,i) longs(j,i)],[lats(j,i) lats(j+1,i)]);
            end
        end
end
