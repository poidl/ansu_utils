function [sns,tns,pns] = depth_ntp_iter(s0,t0,p0,s,t,p)

warning('no check of input dimensions')

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
    F=cast-bottle; % rho-(rho_s+rho'); find corrected surface by finding roots of this term
   
    F_n = F<0;
    
    F_neg=F;
    F_neg(~F_n)=nan;
          
    [min_F, lminr]=min(abs(F_neg)); 
    final=min_F<=delta; % These are points with sufficiently small and negative F. Note that this includes points which are not adjacent to a zero crossing.
            
    zc_F_stable= F_n & circshift(~F_n,-1); % stable zero crossing (F<0 on point above, F>0 on point below); possible candidates for iterative zooming
    zc_F_stable(end,:)=false; % zooming not possible at bottom
    k_zc=nan(1,size(zc_F_stable,2)); % initialize vertical index of zero crossing
    for ii=1:size(zc_F_stable,2); % check if there are multiple stable zero crossings in one vertical cast and remove all except for the shallowest one
        if sum(zc_F_stable(:,ii)>1);
            k_zc(ii)=find(zc_F_stable(:,ii),1,'first'); % remember the vertical index of zero crossing
            zc_F_stable(k_zc(ii)+1:end)=false;
        end
    end
    
    zc_F_stable=any(zc_F_stable,1);
        
    cond1=min_F>delta;
    fr= zc_F_stable & cond1; %  at these horizontal locations we have to increase the vertical resolution before finding the root
    
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


