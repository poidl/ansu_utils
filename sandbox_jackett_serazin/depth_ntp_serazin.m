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
   error('depth_ntp:  Requires two outputs')
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
if(ms*mt*mp == 1 && ns*nt*np ~=1)
    if(ms*mt*mp ~= 1 && ns*nt*np ==1)
        error('depth_ntp: Inputs array dimensions arguments do not agree')
    else
        error('depth_ntp: There must be at least 2 bottles')
    end
end


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
    
    pc0 = p(c) - e*(p(c_d) - p(c))/(e_d - e);
    %Testing undercropping
    if isnan(pc0)
        SAns = NaN;
        CTns = NaN;
        pns=NaN;
        success =1;
    end
    
    %3)Developping some Newton-Raphson iteration
    while success == 0
        iter = iter + 1;
        [SAc0,CTc0] = stp_interp([SA(c),SA(c_d)],[CT(c),CT(c_d)],[p(c),p(c_d)],pc0);
        [sigl,sigu] = sig_vals(SA0,CT0,p0,SAc0,CTc0,pc0);
        ec0 = sigu - sigl;
        p1 = 0.5*(p(c) + pc0);
        ez1 = (e- ec0)/(pc0 - p(c));
        p2 = 0.5*(pc0 + p(c_d));
        ez2 = (ec0 - e_d)/(p(c_d) - pc0);
        r = (pc0 - p1)/(p2 - p1);
        ecz_0 = ez1 + r*(ez2 - ez1);
        if iter == 1
            ecz0 = ecz_0;
        else
            ecz0 = -(ec0 - ec_0)/(pc0 - pc_0);
            if ecz0 == 0
                ecz0 = ecz_0;
            end
        end
        pc1 = pc0 + ec0/ecz0;
        eps = abs(pc1 - pc0);
        %Testing the accuracy
        if abs(ec0) <= 5e-5 && eps <= 5e-3
            SAns = SAc0;
            CTns = CTc0;
            pns = pc0;
            success = 1;
            niter = iter;
        elseif iter > 10
            [SAns,CTns,pns,niter] = e_solve(SA,CT,p,[e e_d],[c c_d],SA0,CT0,p0);
            success = 1;
        else
            pc_0 = pc0;
            ec_0 = ec0;
            pc0 = pc1;
            success = 0;
        end
    end
    
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
    
    %3) Developping some Newton-Raphson iterations
    while success == 0
        iter = iter + 1;
        [SAc0,CTc0] = stp_interp([SA(c_s),SA(c)],[CT(c_s),CT(c)],[p(c_s),p(c)],pc0);
        [sigl,sigu] = sig_vals(SA0,CT0,p0,SAc0,CTc0,pc0);
        ec0 = sigu - sigl;
        p1 = 0.5*(p(c_s) + pc0);
        ez1 = (e_s- ec0)/(pc0 - p(c_s));
        p2 = 0.5*(pc0 + p(c));
        ez2 = (ec0 - e)/(p(c) - pc0);
        r = (pc0 - p1)/(p2 - p1);
        ecz_0 = ez1 + r*(ez2 - ez1);
        if iter == 1
            ecz0 = ecz_0;
        else
            ecz0 = -(ec0 - ec_0)/(pc0 - pc_0);
            if ecz0 == 0
                ecz0 = ecz_0;
            end
        end
        pc1 = pc0 + ec0/ecz0;
        eps = abs(pc1 - pc0);
        %Testing the accuracy
        if abs(ec0) <= 5e-5 && eps <= 5e-3
            SAns = SAc0;
            CTns = CTc0;
            pns = pc0;
            success = 1;
            niter = iter;
        elseif iter > 10
            [SAns,CTns,pns,niter] = e_solve(SA,CT,p,[e_s e],[c_s c],SA0,CT0,p0);
            success = 1;
        else
            pc_0 = pc0;
            ec_0 = ec0;
            pc0 = pc1;
            success = 0;
        end
    end
else
    SAns = NaN;
    CTns = NaN;
    pns=NaN;
end

return