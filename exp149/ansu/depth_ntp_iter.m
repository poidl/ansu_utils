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
   
    [final,fr,k_zc]=root_core(F,delta,stack);
    
    k_zc_3d=k_zc+stack*[0:size(F,2)-1]; % indices of flattened 3d-array where root has been found
    
    sns(inds(final))=s(k_zc_3d(final)); % adjust surface where root has already been found
    tns(inds(final)) =t(k_zc_3d(final));
    pns(inds(final)) =p(k_zc_3d(final));
    inds=inds(fr); % points where surface has not been corrected
    
    if all(~fr) % break out of loop if all roots have been found
        break
    end
    
    k=k_zc_3d(fr);  % indices of flattened 3d-array where vertical resolution must be increased
    
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


