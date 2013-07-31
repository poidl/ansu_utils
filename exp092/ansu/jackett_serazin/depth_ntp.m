function [SAns,CTns,pns] = depth_ntp(SA0,CT0,p0,SA,CT,p)

% depth_ntp                 Absolute Salinity, Conservative Temperature and
%                           in situ Pressure, of the neutral tangent plane
%                           on the neighbouring cast                        
%==========================================================================
% 
% USAGE:  
%  [SAns,CTns,pns] = depth_ntp_v2(SA0,CT0,p0,SA,CT,p)
%
% DESCRIPTION:
%  Find the position in Absolute Salinity, Conservative Temperature and in 
%  situ Pressure where the neutral tangent plane passing through a bottle 
%  intersects a neighbouring cast
%
%  The principle is based on successive approximations
%  - 1: One looks for the simple approximation of neutral tangent plane by
%       finding the right pressure in the neighbouring cast by minizing the  
%       difference between the pressure of the bottle and the cast. It will 
%       be the starting point for the next part 
%                 
%  - 2: One studies then the difference in potential density between the 
%       bottle and the point of the cast. According to the sign, ones looks
%       for the next point denser or less dense. Ones finds in this way an
%       area between a point denser and another less dense for evaluating
%       the position of neutral tangent plane
% 
%  - 3: The position between these two points of the cast is approximate by
%       a Newton-Raphson method.
%
%  INPUT :        
%   SA0  =  the bottle Absolute Salinity                        [ g kg^-1 ]
%   CT0  =  the bottle Conservative Temperature                   [ deg C ]
%   p0   =  the bottle In situ pressure                            [ dbar ]
%   SA   =  vector of cast Absolute Salinities                  [ g kg^-1 ]
%   CT   =  vector of cast Conservative Temperatures              [ deg C ]
%   p    =  vector of cast In situ pressures                       [ dbar ]
% 
%   SA0, CT0 & p0 need to be scalars (dimension 1x1)
%   SA, CT & p need to be vector the dimensions may be Nx1 or 1xN with at 
%   least N>1 
%
%  OUTPUT :       
%   SAns  =  Absolute Salinity of the ntp intersection          [ g kg^-1 ]
%   CTns  =  Conservative Temperature of the ntp intersection     [ deg C ]
%   pns   =  In situ pressure of the ntp intersection              [ dbar ]
%
%  AUTHOR:          
%   David Jacket
%   Modified by Guillaume Serazin 
% 
% VERSION NUMBER: 2.0
%==========================================================================


%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------
if ~(nargin == 6)
   error('depth_ntp:  Requires three inputs')
end
if ~(nargout == 3)
   error('depth_ntp:  Requires three outputs')
end 

[msb,nsb] = size(SA0);
[mtb,ntb] = size(CT0);
[mpb,npb] = size(p0);
[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if(msb*nsb*mtb*ntb*mpb*npb ~= 1)
    error('depth_ntp: Inputs array dimensions arguments do not agree')
end

if (mt ~= ms || mt ~= mp || ns ~= nt || ns ~= np)
    error('depth_ntp: SA and CT must have same dimensions')
end
% if(ms*mt*mp == 1 && ns*nt*np ~=1)
%     if(ms*mt*mp ~= 1 && ns*nt*np ==1)
%         error('depth_ntp: Inputs array dimensions arguments do not agree')
%     else
%         error('depth_ntp: There must be at least 2 bottles')
%     end
% end


%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

n = length(SA);

%1)Looking for the closest pressure to find the starting point
%-------------------------------------------------------------
[~,c]=min(abs(p-p0)); %c for cast


%Evaluating the difference in potential density
[sigl,sigu] = sig_vals(SA0,CT0,p0,SA(c),CT(c),p(c));
e = sigu - sigl;

%Testing exact crossing
%----------------------
if e == 0
    SAns = SA(c);
    CTns = CT(c);
    pns = p(c);
    
%Testing materiality of points
%-----------------------------
elseif isnan(e)
    SAns = NaN;
    CTns = NaN;
    pns=NaN;
    
%Case when starting point less dense than the bottle
%---------------------------------------------------
elseif (e<0&&c<n)
    %Initializing variables
    c_d=c+1;                    %design the next cast deep
    iter = 0;
    success = 0;
    [sigl_d,sigu_d] = sig_vals(SA0,CT0,p0,SA(c_d),CT(c_d),p(c_d));
    e_d=sigu_d - sigl_d;
    
    %2) Looking for the right area
    %While the next point is less dense than the bottle
    %-> going deep
    while(e_d<0&&c_d<n)
        %Testing exact crossing
        if e_d == 0
            SAns = SA(c);
            CTns = CT(c);
            pns = p(c);
            success =1;
            break
        end
        %Reaching the next area deep
        e=e_d;
        c=c_d;
        c_d=c_d+1;
        [sigl_d,sigu_d] = sig_vals(SA0,CT0,p0,SA(c_d),CT(c_d),p(c_d));
        e_d=sigu_d - sigl_d;
    end
    
    % pc0: pressure at which vertical linear interpolation on cast yields bottle
    % dens. refercenced to mid-pressure between bottle and
    % upper-cast-bottle
    pc0 = p(c) - e*(p(c_d) - p(c))/(e_d - e); 
    %Testing undercropping
    if isnan(pc0)
        SAns = NaN;
        CTns = NaN;
        pns=NaN;
        success =1;
    end
    
    %3) find zero crossing with fzero

    su=SA(c); sl=SA(c+1); % linear interpolation of s and ct
    ctu=CT(c); ctl=CT(c+1);
    pu=p(c); pl=p(c+1);
    
    sp=@(p) su+(sl-su)*(p-pu)/(pl-pu); % linear interpolation of s and ct
    ctp=@(p) ctu+(ctl-ctu)*(p-pu)/(pl-pu);
    
    pc=@(pmid) 2*pmid-p0; % pressure on cast, as a function of mid-pressure (and the paramter p0; bottle pressure)
    fac=1./sqrt(gsw_rho(su,ctu,pu)^2+gsw_rho(sl,ctl,pl)^2); % scale to avoid having to set tolerance ?
    %fac=1;
    func_normalized=@(pmid)  fac.*(gsw_rho(sp(pc(pmid)), ctp(pc(pmid)), pmid) -gsw_rho(SA0,CT0,pmid));
    
    proot=fzero(func_normalized, 0.5*(p0+[pu,pl]));
    
    pns=pc(proot);
    SAns= su+(sl-su)*(pns-pu)/(pl-pu);
    CTns=ctu+(ctl-ctu)*(pns-pu)/(pl-pu);
    
    
%Case when starting point is denser than the bottle
%--------------------------------------------------
elseif (e>0&&c>1)
    c_s=c-1;    %design the next cast shallow
    success = 0;
    iter = 0;
    [sigl_u,sigu_u] = sig_vals(SA0,CT0,p0,SA(c_s),CT(c_s),p(c_s));
    e_s=sigu_u - sigl_u;
    
    %2) Looking for the right area
    %While the next point is denser than the bottle :
    %-> going shallow
    while(e_s>0&&c_s>1)
        %Testing exact crossing
        if e_s == 0
            SAns = SA(c);
            CTns = CT(c);
            pns = p(c);
            success =1;
            break
        end
        %Reaching the next point shallow
        e=e_s;
        c=c_s;
        c_s=c_s-1;  
        [sigl_u,sigu_u] = sig_vals(SA0,CT0,p0,SA(c_s),CT(c_s),p(c_s));
        e_s=sigu_u - sigl_u;
    end
    
    pc0 = p(c_s) - e_s*(p(c) - p(c_s))/(e - e_s);
    %Testing outcropping
    if pc0<-1.5
        SAns = NaN;
        CTns = NaN;
        pns=NaN;
        success =1;
    end
    
    %3) find zero crossing with fzero

    su=SA(c-1); sl=SA(c); % linear interpolation of s and ct
    ctu=CT(c-1); ctl=CT(c);
    pu=p(c-1); pl=p(c);
    
    sp=@(p) su+(sl-su)*(p-pu)/(pl-pu); % linear interpolation of s and ct
    ctp=@(p) ctu+(ctl-ctu)*(p-pu)/(pl-pu);
    
    pc=@(pmid) 2*pmid-p0; % pressure on cast, as a function of mid-pressure (and the paramter p0; bottle pressure)
    fac=1./sqrt(gsw_rho(su,ctu,pu)^2+gsw_rho(sl,ctl,pl)^2); % scale to avoid having to set tolerance ?
    %fac=1;
    func_normalized=@(pmid)  fac.*(gsw_rho(sp(pc(pmid)), ctp(pc(pmid)), pmid) -gsw_rho(SA0,CT0,pmid));
    
    proot=fzero(func_normalized, 0.5*(p0+[pu,pl]));
    
    pns=pc(proot);
    SAns= su+(sl-su)*(pns-pu)/(pl-pu);
    CTns=ctu+(ctl-ctu)*(pns-pu)/(pl-pu);
    
else
    SAns = NaN;
    CTns = NaN;
    pns=NaN;
end

return

