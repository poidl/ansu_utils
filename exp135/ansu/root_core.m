function [final,fr,k_zc]=root_core(F,delta,stack)

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

end
