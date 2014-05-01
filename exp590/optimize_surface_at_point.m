function [sns,ctns,pns] = optimize_surface_at_point(s,ct,p,point)

% construct an optimized approximately neutral surface intersecting with 
% point=[p0 ilat ilon], whose horizontal coordinates coincide with a grid point of
% the hydrography, and whose pressure is arbitrary

% ilon:  zonal index with respect to hydrography grid
% ilat:  meridional index with respect to hydrography grid
% p0: pressure at point

% The output surface is one single connected surface. 

[nz,ny,nx]=size(s);
p0=point(1);
ilat=point(2);
ilon=point(3);

pns=p0+zeros(ny,nx);
sns=var_on_surf_stef(s,p,pns);
ctns=var_on_surf_stef(ct,p,pns);


pns(isnan(sns))=nan;

% Although the output surface is one single connected surface, the input to
% the optimization is a global isobar, consisting of possibly disconnected 
% surfaces.
% If only a single connected surface is supplied, this can lead to a lot
% of "append and optimize" iterations in the automated lateral extension of
% the surface. Consider the following example: We want to know the
% optimized surface intersection with the point=[1000 dbar at 16 South, 188
% East] and start with the single connected 1000 dbar isobar intersecting
% with the point. This surface may be disconnected from the artic ocean (if it
% is deeper than Denmark Strait etc.), depending on the horizontal
% resolution). During adjustment, the surface will strongly shoal at high
% latitudes and will be automatically extended through Denmark Strait etc.
% However, only one grid point at a time is appended between consecutive
% iterations, and it will take many iterations to fill out the entire
% Arctic Ocean. 
% This can be avoided by supplying a global isobar to the
% optimization. The Arctic will be adjusted as an isolated Basin in the
% first iteration, and then automatically joined to the Atlantic.
% TODO: This approach may be disadvantagous for surfaces that do not "grow"
% considerably.

[sns,ctns,pns] = optimize_surface_exact(s,ct,p,sns,ctns,pns);

keyboard
istation=ilat+ny*(ilon-1);

setnan=true(size(sns));
regions=find_regions(sns);
for iregion=1:length(regions)
    region=regions{iregion};
    if ismember(istation,region)
        setnan(region)=false;        
        pns(setnan)=nan;
        sns(setnan)=nan;
        ctns(setnan)=nan;      
    end
end 

end
