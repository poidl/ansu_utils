function [rkx,rky]=delta_rhokappa(s,ct,p)
omega_user_input;

% forward differences
rk=gsw_rho(s,ct,p).*gsw_kappa(s,ct,p); % rk = rho*kappa = drho/dp

rkx=circshift(rk, [0 0 -1])-rk;
rky=circshift(rk, [0 -1 0])-rk;

if ~zonally_periodic
    rkx(:,:,end)=nan;
end
rky(:,end,:)=nan;


% on eastern (northern) boundaries, make a sloppy backward difference (the
% objective is to 'lose' as few grid points as possible, even it makes the 
% fd scheme inconsistent)

inx=~isnan(rk) & isnan(rkx);
iny=~isnan(rk) & isnan(rky);

rkx_=rk-circshift(rk, [0 0 1]);
rky_=rk-circshift(rk, [0 1 0]);
if ~zonally_periodic
    rkx_(:,:,1)=nan;
end
rky_(:,end,:)=nan;

rkx(inx)=rkx_(inx);
rky(iny)=rky_(iny);

inx2=~isnan(rk) & isnan(rkx);
iny2=~isnan(rk) & isnan(rky);

% save_netcdf(rk,'rk','data/rk.nc');
% save_netcdf(rkx,'rkx','data/rkx.nc');
% save_netcdf(rky,'rky','data/rky.nc');
% keyboard
end