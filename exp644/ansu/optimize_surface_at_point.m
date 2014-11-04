function [sns,ctns,pns] = optimize_surface_at_point(s,ct,p,lon,lat,point)

    % construct an optimized approximately neutral surface intersecting with 
    % point=[p0 ilat ilon], whose horizontal coordinates coincide with a grid point of
    % the hydrography, and whose pressure is arbitrary

    % ilon:  zonal index with respect to hydrography grid
    % ilat:  meridional index with respect to hydrography grid
    % p0: pressure at point

    % The output surface is one single connected surface. 

    omega_user_input;

    [nz,ny,nx]=size(s);
    p0=point(1);
    ibb=point(2);

    % as initial surface, we choose an iso-surface of potential density with reference pressure p0 
    sig=gsw_rho(s,ct,p0+zeros(size(s))); 
    sigc=sig(:,ibb);
    pc=p(:,ibb);
    si=var_on_surf_stef(sigc,pc,p0); % pot. dens. value of initial surface
    si=si+zeros(ny,nx);

    sns=var_on_surf_stef(s,sig,si);
    ctns=var_on_surf_stef(ct,sig,si);
    pns=var_on_surf_stef(p,sig,si);
    
    if abs(pns(ibb)-p0)>1e-8
        error('problem')
    end
    
%     % start from isobar
%     pns= p0+zeros(ny,nx);
%     sns=var_on_surf_stef(s,p,pns);
%     ctns=var_on_surf_stef(ct,p,pns);
%     pns(isnan(sns))=nan;

    [dx,dy,dz]=get_dx(lon,lat,p);
    dx=squeeze(dx(1,:,:));
    dy=squeeze(dy(1,:,:));
    save([datapath,'dxdy.mat'], 'dx', 'dy') 
    save([datapath,'ibb.mat'], 'ibb')

    if no_land_mask
        save([datapath,'latlon.mat'],'lat','lon')
    end

    % we only keep the one single (connected) surface on which p0 is located.
    % Attention if no_land_mask=true !
    [sns,ctns,pns] = get_connected(sns,ctns,pns,ibb);
    
    [sns,ctns,pns,errx,erry] = optimize_surface_exact(s,ct,p,sns,ctns,pns);

end

