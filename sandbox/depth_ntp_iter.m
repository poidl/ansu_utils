function [sns,tns,pns] = depth_ntp_iter(s0,t0,p0,s,t,p)

s=s(:,:);
t=t(:,:);
p=p(:,:);

%ii=bsxfun(@times,1:size(s,2),ones(size(s,1),1));
% 
% s0=s0(ii);
% t0=t0(ii);
% p0=p0(ii);

%pmid=0.5*(p0+p);

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
        ii=bsxfun(@times,1:size(s,2),ones(stack,1));
        s0_stacked=s0(ii);
        t0_stacked=t0(ii);
        p0_stacked=p0(ii);
    end
    
    inds=inds(fr); % points where surface has not been corrected
    s0_stacked=s0_stacked(:,fr);
    t0_stacked=t0_stacked(:,fr);
    p0_stacked=p0_stacked(:,fr);
    
    pmid=0.5*(p0_stacked+p);
    bottle=gsw_rho(s0_stacked,t0_stacked,pmid);

    cast=gsw_rho(s(:,:),t(:,:),pmid); % 3-d density referenced to pressure of the current surface
    dr=cast-bottle; % rho-(rho_s+rho'); find corrected surface by finding roots of this term
    
    dr_p = dr>0;
    dr_n = dr<0;
    
    zc = any(dr_p,1) & any(dr_n,1); % horizontal indices of locations where zero-crossing occurs
    [min_dr, lminr]=min(abs(dr));
    cond1=min_dr>delta;
    final=min_dr<=delta; % at these horizontal indices root has been found
    fr= zc & cond1; %  at these horizontal locations we have to increase the vertical resolution before finding the root
    
    lminr=lminr+stack*[0:size(dr_n,2)-1];
    lminr=lminr(final);
    
    sns(inds(final))=s(lminr); % adjust surface where root has already been found
    tns(inds(final)) =t(lminr);
    pns(inds(final)) =p(lminr);
    
    
    if all(~fr) % break out of loop if all roots have been found
        break
    end
    
    k=sum(dr_n,1); % find indices of flattened 3d-array where vertical resolution must be increased
    k=k+stack*[0:size(dr_n,2)-1];
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


