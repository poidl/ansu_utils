function [SAns,CTns,pns,iter] = e_solve(SA,CT,p,e,k,SA0,CT0,p0)

%	Find the zero of the e function using a bisection method
%
%  Usage :         [SAns,CTns,pns,iter] = e_solve(SA,CT,p,e,k,SA0,CT0,p0)
%
%	Input :		SA		array of cast Absolute Salinities
%				CT		array of cast Conservative Temperatures
%				p		array of cast pressures
%				e		array of cast e values
%				k		interval (k-1, k) contains the zero
%				SA0		the bottle Absolute Salinity
%				CT0		the bottle Conservative Temperature
%				p0		the bottle pressure
%
%	Output :	SAns		Absolute Salinity of e zero
%				CTns		Conservative Temperature of e zero
%				pns		    pressure of e zero
%
%	Units :	    Absolute Salinities	        g/kg (TEOS-10)
%				Conservative Temperatures	degrees C (ITS-90)
%				pressures	dbar
%  DRJ on 17/06/03

kl=k(1);
ku=k(2);
pl = p(kl);
el = e(1);
pu = p(ku);
eu = e(2);

iter = 0;
success = 0;

while success == 0
    iter = iter + 1;
    pm = 0.5*(pl + pu);
    [SAm, CTm] = stp_interp([SA(kl),SA(ku)],[CT(kl),CT(ku)],[p(kl),p(ku)],pm);
    [sigl, sigu] = sig_vals(SA0,CT0,p0,SAm,CTm,pm);
    em = sigu - sigl;
    if el*em < 0
        pu = pm;
        eu = em;
    elseif em*eu < 0
        pl = pm;
        el = em;
    elseif em == 0
        SAns = SAm;
        CTns = CTm;
        pns = pm;
        success = 1;
    end
    if success == 0
        if (abs(em) <= 5e-5) && (abs(pu-pl) <= 5e-3)
            SAns = SAm;
            CTns = CTm;
            pns = pm;
            success = 1;
        elseif iter <= 20
            success = 0;
        else
            disp('WARNING in e-solve')
            disp(['iter: ', int2str(iter), '  em: ', num2str(em), '  dp: ', num2str(abs(pu-pl))])
            SAns = NaN;
            CTns = NaN;
            pns = NaN;
            success = 0;
            return
        end
    end
end


return