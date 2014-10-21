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

    % This is on crack; sig may be non-monotonous and iso-surfaces may be
    % branched. Better start from isobar.

    % as initial surface, we choose an iso-surface of potential density with reference pressure p0 
    sig=gsw_rho(s,ct,p0+zeros(size(s))); 
    sigc=sig(:,ilat,ilon);
    pc=p(:,ilat,ilon);
    si=var_on_surf_stef(sigc,pc,p0); % pot. dens. value of initial surface
    si=si+zeros(ny,nx);

    sns=var_on_surf_stef(s,sig,si);
    ctns=var_on_surf_stef(ct,sig,si);
    pns=var_on_surf_stef(p,sig,si);
    
    if abs(pns(ilat,ilon)-p0)>1e-8
        error('problem')
    end
    
%     % start from isobar
%     pns= p0+zeros(ny,nx);
%     sns=var_on_surf_stef(s,p,pns);
%     ctns=var_on_surf_stef(ct,p,pns);
%     pns(isnan(sns))=nan;
    
    istation=ilat+ny*(ilon-1);

    % we only keep the one single (connected) surface on which p0 is located.
    % Attention if no_land_mask=true !
    [sns,ctns,pns] = get_connected(sns,ctns,pns,istation);
    
    [sns,ctns,pns,errx,erry] = optimize_surface_exact(s,ct,p,sns,ctns,pns);

end

