function [sns,tns,pns] = depth_ntp_iter(s0,t0,p0,s,t,p)

%warning('no check of input dimensions')

s=s(:,:);
t=t(:,:);
p=p(:,:);

zi=size(s,1);
yixi=size(s,2);
refine_ints=100;

inds=1:yixi;
fr=true(1,yixi);
delta=1e-9;

pns = nan(1,yixi);
sns = nan(1,yixi);
tns = nan(1,yixi);

cnt=0;
while 1
    cnt=cnt+1;
    
    if cnt==1 % in first iteration pns_l is stacked vertically zi times, after that it is stacked refine_ints times
        stack=zi;
    elseif cnt==2
        stack=refine_ints+1;
    end
    if cnt==1 | cnt==2
        ii=bsxfun(@times,1:yixi,ones(stack,1));
        s0_stacked=s0(ii);
        t0_stacked=t0(ii);
        p0_stacked=p0(ii);
    end
    
    s0_stacked=s0_stacked(:,fr);
    t0_stacked=t0_stacked(:,fr);
    p0_stacked=p0_stacked(:,fr);
    
    pmid=0.5*(p0_stacked+p);
    bottle=gsw_rho(s0_stacked,t0_stacked,pmid);

    cast=gsw_rho(s(:,:),t(:,:),pmid); % 3-d density referenced to pmid
    F=cast-bottle; 
   
    F_p = F>0;
    F_n = F<0;
    
    % TODO: a double-zero-crossing could arise due to linear interpolation for
    % values that are close to 0 and of equal sign in F?
    
    zc_F_stable= F_n & circshift(F_p,-1); % stable zero crossing (F<0 at point and F>0 on point below);
    zc_F_stable(end,:)=false; % discard bottom (TODO: should check if bottom point is negative, sufficiently close to zero and has a negative point above it)
    
    k_zc=nan(1,size(zc_F_stable,2)); % initialize vertical index of zero crossing
    any_zc_F_stable=any(zc_F_stable,1);
    for ii=find(any_zc_F_stable); % check if there are multiple stable zero crossings
        zc1=find(zc_F_stable(:,ii),1,'first'); % get the vertical index of shallowest zero crossing
        zc_F_stable(zc1+1:end,ii)=false; % remove additional crossings
        k_zc(ii)=zc1; % remember the vertical index of shallowest zero crossing
    end
    
    F_neg=F;
    F_neg(~zc_F_stable)=nan; %  Matrix of negative F values which lie above zero crossings.
          
    [min_F, lminr]=min(abs(F_neg)); 
    final=(min_F<=delta); % These are points with sufficiently small F.
    
    cond1=min_F>delta;
    fr= any_zc_F_stable & cond1; %  at these horizontal locations we have to increase the vertical resolution before finding the root
    
    lminr=lminr+stack*[0:size(F_n,2)-1];
    lminr=lminr(final); 
    
    sns(inds(final))=s(lminr); % adjust surface where root has already been found
    tns(inds(final)) =t(lminr);
    pns(inds(final)) =p(lminr);
    inds=inds(fr); % points where surface has not been corrected
    
    if all(~fr) % break out of loop if all roots have been found
        break
    end
    
    k=k_zc; % find indices of flattened 3d-array where vertical resolution must be increased
    k=k+stack*[0:size(F_n,2)-1];
    k=k(fr);
    
    ds_ =  ( s(k+1) - s(k))/refine_ints; % increase resolution in the vertical
    dt_ = (t(k+1) - t(k))/refine_ints;
    dp_ =  (p(k+1) - p(k))/refine_ints;
    
    ds_ =bsxfun(@times, ds_, [0:refine_ints]');
    dt_ = bsxfun(@times, dt_, [0:refine_ints]');
    dp_ = bsxfun(@times, dp_, [0:refine_ints]');
    
    s =  bsxfun(@plus,s(k),ds_);
    t =  bsxfun(@plus,t(k),dt_);
    p =  bsxfun(@plus,p(k),dp_);
    
end


