function [erx,ery]=no_land_mask_disconnect_rule(erx,ery)
% in case there is no land mask (e.g. in a climatology data set with too
% coarse resolution to represent land), set erx and ery to nan at the
% appropriate points. The rule where to set nans is defined here. 

% Example below refers to the Jackett/McDougall 96 data set, and disconnects North Indian from North Pacific
% (Thailand/Malaysia/Indonesia not resolved) and North Pacific from North Atlantic (Central America not resolved)

%%%%% Example start
% omega_user_input;
% load([datapath,'no_land_mask.mat']) % load output of gamma_ocean_and_n(s,ct,p,lon,lat);
% [ny,nx]=size(erx);
% test=ocean(:).*circshift(ocean(:),-ny);
% en_in_other_basin= (test==5); % eastern neighbour in other basin
% erx(en_in_other_basin)=nan;
% test=ocean(:).*circshift(ocean(:),-1);
% nn_in_other_basin= (test==5); % northern neighbour in other basin
% ery(nn_in_other_basin)=nan;
%%%%% Example end

omega_user_input;

load([datapath,'no_land_mask.mat'])

[ny,nx]=size(erx);

test=ocean(:).*circshift(ocean(:),-ny);
en_in_other_basin= (test==5); % eastern neighbour in other basin
erx(en_in_other_basin)=nan;

test=ocean(:).*circshift(ocean(:),-1);
nn_in_other_basin= (test==5); % northern neighbour in other basin
ery(nn_in_other_basin)=nan;
        

