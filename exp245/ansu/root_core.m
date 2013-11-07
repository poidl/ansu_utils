function [s,ct,p,sns_out,ctns_out,pns_out,inds, fr, dobreak]=root_core(F,delta,stack,inds,refine_ints,s,ct,p,sns_out,ctns_out,pns_out);
    dobreak=false;
% final: boolean array of horizontal positions, indicating that root has been found
% fr: boolean array of horizontal positions, indicating that there is a stable zero crossing and root has not been found yet (zoom here)
% k_zc: a vector containing the vertical position of the zero crossing for each horizontal position

    F_p = F>=0;
    F_n = F<0;
    
    % TODO: a double-zero-crossing could arise due to linear interpolation for
    % values that are close to 0 and of equal sign in F?
    
    zc_F_stable= F_n & circshift(F_p,-1); % stable zero crossing (F<0 at point and F>0 on point below);
    zc_F_stable(end,:)=false; % discard bottom (TODO: should check if bottom point is negative, sufficiently close to zero and has a negative point above it)
    
    cs=cumsum(zc_F_stable,1)+1;
    cs(cs~=1)=0;
    k_zc=sum(cs,1)+1;% vertical index of shallowest stable zero crossing
    any_zc_F_stable=any(zc_F_stable,1);
    k_zc(~any_zc_F_stable)=1; % dummy to avoid zeros as indices

    F_neg=F(k_zc+stack*[0:size(F,2)-1]); % value of F above the shallowest stable zero crossing (or dummy if there is no stable zero crossing)
    F_neg(~any_zc_F_stable)=nan; % remove dummy indices
    
    final=(abs(F_neg)<=delta); % These are points with sufficiently small F.

    cond1=abs(F_neg)>delta;
    fr= any_zc_F_stable & cond1; %  at these horizontal locations we have to increase the vertical resolution before finding the root    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    k_zc_3d=k_zc+stack*[0:size(F,2)-1]; % indices of flattened 3d-array where root has been found   
    
    sns_out(inds(final))=s(k_zc_3d(final)); % adjust surface where root has already been found
    ctns_out(inds(final)) =ct(k_zc_3d(final));
    pns_out(inds(final)) =p(k_zc_3d(final));
    inds=inds(fr); % points where surface has not been corrected
    
    if all(~fr) % break out of loop if all roots have been found
        dobreak=true;
        return
    end
    
    k=k_zc_3d(fr);  % indices of flattened 3d-array where vertical resolution must be increased
    
    ds_ =  ( s(k+1) - s(k))/refine_ints; % increase resolution in the vertical
    dct_ = (ct(k+1) - ct(k))/refine_ints;
    dp_ =  (p(k+1) - p(k))/refine_ints;
    
    ds_ =bsxfun(@times, ds_, [0:refine_ints]');
    dct_ = bsxfun(@times, dct_, [0:refine_ints]');
    dp_ = bsxfun(@times, dp_, [0:refine_ints]');
    
    s =  bsxfun(@plus,s(k),ds_);
    ct =  bsxfun(@plus,ct(k),dct_);
    p =  bsxfun(@plus,p(k),dp_);
    
end
