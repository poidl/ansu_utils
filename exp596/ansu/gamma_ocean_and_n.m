function [ocean, n] = gamma_ocean_and_n(SA,CT,p,long,lat)

% gamma_ocean_and_n
%==========================================================================
%
% USAGE:  
%
%==========================================================================

% the ocean array
dim_SA = size(SA);
dim_lat = size(lat);
dim_long = size(long);

if length(dim_SA) == 3
    if length(dim_lat) == 1 & length(dim_long) == 1
        ocean = gamma_ocean_3d(long,lat,squeeze(SA(1,:,:)));
    elseif length(dim_lat) == 1 & length(dim_long) == 2
        ocean = gamma_ocean_3d(squeeze(long(1,:)),squeeze(lat(:,1)),squeeze(SA(1,:,:)));
    elseif length(dim_lat) == 1 & length(dim_long) == 1
        ocean = gamma_ocean_3d(squeeze(long(1,1,:)),squeeze(lat(1,1,:)),squeeze(SA(1,:,:)));
    else
        ocean = gamma_ocean_3d(unique(sort(long(:))),unique(sort(lat(:))),squeeze(SA(1,:,:)));
    end
    
elseif length(dim_SA) == 2 && dim_SA(2) > 1
    nx = dim_SA(2);
    ocean = nan(1,nx);
    for i = 1:nx
        ocean(i) = gamma_ocean_2d(long(i),lat(i));
    end
else
    error('****    not a section or ocean of data    ****')
end

% the n array
n = ones(size(SA));
data = SA.*CT.*p;
Inan = find(isnan(data));
if ~isempty(Inan)
    n(Inan) = nan;
end
n = reshape(nansum(n),size(ocean));

end

%##########################################################################

function ocean = gamma_ocean_3d(long,lat,z)

% gamma_ocean_3d
%==========================================================================
%
% USAGE:  
%  ocean = gamma_ocean_3d(longs,lats,z)
%
% DESCRIPTION:
%  gamma_ocean_3d      ocean of longs/lats observations
%
% INPUT: 
%  longs - vector (nx) of longitudes [0,360]
%  lats  - vector (ny) of latitudes [-90,90]
%  z     - array (ny,nx) describing the ocean extent
%
% OUTPUT:   
%  ocean - array (ny,nx) of ocean values
%           0    land,
%           1    North Pacific,
%           2    South Pacific
%           4    South Indian
%           5    North Atlantic, Hudson Bay, Mediterranean Sea, North Indian
%                Red Sea, persian gulf, Baltic sea, black/caspian seas
%           6    South Atlantic
%           7    Atlantic / Artic
%           8    Southern Equatorial Pacific, 
%           9    
%           10   
% Author:    	
%  David Jackett 7/7/98
%
%==========================================================================

nx = length(long);
ny = length(lat);
ocean = nan(ny,nx);

for j = 1:ny
    for i = 1:nx
        if ~isnan(z(j,i))
            ocean(j,i) = gamma_ocean_2d(long(i),lat(j));
        end
    end
end

end

%##########################################################################

function ocean = gamma_ocean_2d(long,lat)

% gamma_ocean_2d
%==========================================================================
%
% USAGE:  
%  ocean = gamma_ocean0(long,lat)
%
%   gamma_ocean0         ocean of single long/lat observation
%
% Input: 
%  long  - 	longitude [0,360]
%  lat   - 	latitude [-90,90]
%
% Output:    
%  ocean - 	0  land,
%           1-6  main oceans,
%           7    Arctic
%           8    Med
%
% Author:
%  David Jackett  7/7/98
%
%==========================================================================


po_long = [100, 140, 240, 260, 272.59, 276.5, 278.65, 280.73, 295.217, 290, ...
    300, 294, 290, 146, 146, 133.9, 126.94, 123.62, 120.92, 117.42, 114.11, ...
    107.79, 102.57, 102.57, 98.79, 100];

po_lat = [20, 66, 66, 19.55, 13.97, 9.6, 8.1, 9.33, 0, -52, -64.5, -67.5, ...
    -90, -90, -41, -12.48, -8.58, -8.39, -8.7, -8.82, -8.02, -7.04, -3.784, ...
    2.9, 10, 20];

na_long = [260, 272, 285, 310, 341.5, 353.25, 356, 360, 354.5, 360, 360, ...
    295.217, 280.73, 278.65, 276.5, 272.59, 260, 260];

na_lat = [60, 72, 82, 81.75, 65, 62.1, 56.5, 44, 36, 20, 0, 0, 9.33, 8.1, ...
    9.6, 13.97, 19.55, 60];

i_pacific = inpolygon(long,lat,po_long,po_lat);

% pacific ocean
if i_pacific == 1
    if lat <= -15
        ocean = 5;%2;
    elseif lat > -15 & lat <= 0
        ocean = 8;
    else
        ocean = 1;
    end
    
% indian ocean
elseif 20<=long && long<=150 && -90<=lat && lat<=30
    if lat <= 0
        ocean = 4;
    else
        ocean = 5;
    end
    
% atlantic ocean
elseif (0<=long && long<=20 && -90<=lat && lat<=90) || ...
        (20<=long && long<=40 && 28<=lat && lat<=44) || ...
        (260<=long && long<=360 && -90<=lat && lat<=90)
    i_natlantic = inpolygon(long,lat,na_long,na_lat);
%     if lat <= 0
%         ocean = 5;
%     else
    if (i_natlantic==1) | (0<=long && long<=15 && 0<lat && lat<=10) | lat <= 0
        ocean = 5;
    else
        ocean = 7;
    end
    
% arctic ocean
elseif 0<=long && long<=360 && 64<lat && lat<=90
    ocean = 7;
else
    ocean = 7;
end

% % red sea
% if 31.25<=long && long<=43.25 && 13<=lat && lat<=30.1
%     ocean = 5;
%     
% % persian gulf
% elseif 40<=long && long<=56 && 22<=lat && lat<=32
%     ocean = 5;
%     
% % baltic sea
% elseif 12<=long && long<=33.5 && 50<=lat && lat<=60
%     ocean = 9;
% elseif 16<=long && long<=33.5 && 60<=lat && lat<=66
%     ocean = 9;
%     
% % black/caspian seas
% elseif 40<=long && long<=112 && 30<=lat && lat<=60
%     ocean = 10;
% elseif 27.5<=long && long<=40 && 40<=lat && lat<=50
%     ocean = 10;
%     
% % african lakes
% elseif 29<=long && long<=36 && -16<=lat && lat<=4
%     ocean = 10;
%     
% % hudson bay
% elseif 266<=long && long<=284 && 40<=lat && lat<=70
%     ocean = 5;
%    

% mediterranean sea
if 0<=long && long<=45 && 28<=lat && lat<=47
    ocean = 5;
elseif 355<=long && long<=360 && 34<=lat && lat<=42
    ocean = 5;
end

end

%##########################################################################

