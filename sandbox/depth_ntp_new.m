function [sns,tns,pns] = depth_ntp(s0,t0,p0,s,t,p)

%%  Find the position where the neutral tangent plane passing through a bottle 
%%  intersects a neighbouring cast 
%%
%%  Usage :        [sns,tns,pns] = depth_ntp(s0,t0,p0,s,t,p)
%%
%%  Input :        s0    the bottle salinity
%%                 t0    the bottle conservative temperature
%%                 p0    the bottle pressure
%%                 s     vector of cast salinities
%%                 t     vector of cast conservative temperatures
%%                 p     vector of cast pressures
%%
%%  Output :       sns   salinity of the ntp intersection
%%                 tns   conservative temperature of the intersection
%%                 pns   pressure of the intersection
%%
%%  Units :        salinities	  psu (IPSS-78)
%%                 temperatures   degrees C (IPS-90)
%%                 pressures      db

%%  DRJ on 17/06/03


n = length(s); e = zeros(n,1);

%		find the bottle pairs containing a crossing

ncr = 0;
for k = 1:n
  [sigl,sigu] = sig_vals(s0,t0,p0,s(k),t(k),p(k));
  e(k) = sigu-sigl;
  if k>1    
    if e(k-1)==0                   %  an exact crossing at the k-1 bottle
      ncr = ncr+1; sns = s(k-1); tns = t(k-1); pns = p(k-1);
    elseif e(k)*e(k-1)<0           %  a crossing between k-1 and k bottles
      ncr = ncr+1;
                                   %  some Newton-Raphson iterations
      pc0 = p(k-1)-e(k-1)*(p(k)-p(k-1))/(e(k)-e(k-1));
      iter = 0; success = 0;
      if success==0
        su=s(k-1); sl=s(k); % linear interpolation of s and ct
        ctu=t(k-1); ctl=t(k);
        pu=p(k-1); pl=p(k);
        
        sp=@(p) su+(sl-su)*(p-pu)/(pl-pu); % linear interpolation of s and ct
        ctp=@(p) ctu+(ctl-ctu)*(p-pu)/(pl-pu);
        
        pc=@(pmid) 2*pmid-p0; % pressure on cast, as a function of mid-pressure (and the paramter p0; bottle pressure)
        fac=1./sqrt(gsw_rho(su,ctu,pu)^2+gsw_rho(sl,ctl,pl)^2); % scale to avoid having to set tolerance ?
        %fac=1;
        func_normalized=@(pmid)  fac.*(gsw_rho(sp(pc(pmid)), ctp(pc(pmid)), pmid) -gsw_rho(s0,t0,pmid));
        
        proot=fzero(func_normalized, 0.5*(p0+[pu,pl]));
        
        pns=pc(proot);
        sns= su+(sl-su)*(pns-pu)/(pl-pu);
        tns=ctu+(ctl-ctu)*(pns-pu)/(pl-pu);

      end
    end
  end
  if k==n&e(k)==0                  %  the last bottle
    ncr = ncr+1;
	sns = s(k); tns = t(k); pns = p(k);
  end
end

                                   %  multiple and no crossings

if ncr==0
  if e(1)>0                                 %  outcropping
    sns = nan; tns = nan; pns = nan;
    %sns = -99.1; tns = -99.1; pns = -99.1;
  else                                      %  undercropping
    sns = nan; tns = nan; pns = nan;
    %sns = -99.2; tns = -99.2; pns = -99.2;
  end
%elseif ncr>=2                               %  multiple crossings
%	sns = -99.3; tns = -99.3; pns = -99.3;
end

    
return