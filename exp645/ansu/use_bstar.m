

function [drhodx,drhody,regions,b]=use_bstar(drhodx,drhody,pns,s,ct,p)
    omega_user_input;
    %keyboard
    
    g=9.81;
    
    [zi,yi,xi]=size(s);
    [n2,pmid]=gsw_Nsquared(s,ct,p);
    
    n2=reshape(n2,[zi-1,yi,xi]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     pmid=reshape(pmid,[zi-1,yi,xi]);
%  
%     % regrid onto drhox and drhoy grids
%     n2x=  0.5*(circshift(n2,   [0 0 -1])+n2);
%     n2y=  0.5*(circshift(n2,   [0 -1 0])+n2);
%     pmidx=0.5*(circshift(pmid, [0 0 -1])+pmid);
%     pmidy=0.5*(circshift(pmid, [0 -1 0])+pmid);
%  
%     [rkx,rky]=delta_rhokappa(s,ct,p); 
%     
%     rkx=0.5*(circshift(rkx,   [-1 0 0])+rkx); % regrid onto vertical n2 grid
%     rky=0.5*(circshift(rky,   [-1 0 0])+rky);
%     rkx=rkx(1:end-1,:,:);
%     rky=rky(1:end-1,:,:);
% 
%     lnbx=(g^2./n2x).*rkx;
%     lnby=(g^2./n2y).*rky;    
%     
%     pnsx= 0.5*(circshift(pns,    [0 -1])+pns);
%     pnsy= 0.5*(circshift(pns,    [-1 0])+pns);
%     
%     lnbx=var_on_surf_stef(lnbx,pmidx,pnsx);
%     lnby=var_on_surf_stef(lnby,pmidy,pnsy);
% 
%     [regions]=remove_points(lnbx,lnby,pns);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    n2=cat(1,n2(1,:,:),n2); % sloppy
    [rkx,rky]=delta_rhokappa(s,ct,p); % sloppy
    
    lnbx=(g^2./n2).*rkx;
    lnby=(g^2./n2).*rky;
    
    lnbx=var_on_surf_stef(lnbx,p,pns);
    lnby=var_on_surf_stef(lnby,p,pns);
    
    lnbx=0.5*(circshift(lnbx,    [0 -1])+lnbx);
    lnby=0.5*(circshift(lnby,    [-1 0])+lnby);

    %lnbx=over_dA(lnbx,'i');
    %lnby=over_dA(lnby,'j');
    
    [regions]=remove_points(lnbx,lnby,pns);
    lnb=solve_lsqr(regions,lnbx,lnby);
    %lnb=times_dA(lnb);
    
    b=exp(lnb);
    
    bx=0.5*(circshift(b,[0 -1])+b);
    by=0.5*(circshift(b,[-1 0])+b);

    drhodx=bx.*drhodx;
    drhody=by.*drhody;
    
end
