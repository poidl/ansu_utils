function [rkx,rky]=delta_rhokappa(s,ct,p)
user_input;

rk=gsw_rho(s,ct,p).*gsw_kappa(s,ct,p); % rk = rho*kappa = drho/dp

rkx=circshift(rk, [0 0 -1])-rk;
rky=circshift(rk, [0 -1 0])-rk;

if ~zonally_periodic
    rkx(:,:,end)=nan;
end
rky(:,end,:)=nan;

end