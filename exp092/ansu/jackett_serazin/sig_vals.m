function [sig1,sig2] =  sig_vals(SA1,CT1,p1,SA2,CT2,p2)

%	Computes the sigma values of two neighbouring bottles w.r.t. the mid pressure
%
%   Usage :         [sig1,sig2] =  sig_vals(SA1,CT1,p1,SA2,CT2,p2)
%
%	Input :			SA1, SA2         bottle salinities
%					CT1, CT2         bottle Conservative Temperatures
%					p1, p2           bottle pressures
%
%	Output :		sig1, sig2     bottle potential density values
%
%	Units :			SA      g/kg (TEOS-10)
%					CT      degrees C (IPTS-90)
%					pressure      dbar
%                   density       kg m-3
%  DRJ on 17/06/03

pmid = 0.5*(p1 + p2);
sig1 = gsw_rho_CT(SA1,CT1,pmid) - 1000;
sig2 = gsw_rho_CT(SA2,CT2,pmid) - 1000;

return