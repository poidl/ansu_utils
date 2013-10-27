function [sns,ctns,pns] = depth_ntp_iter(s0,ct0,p0,s,ct,p)

%warning('no check of input dimensions')

s=s(:,:);
ct=ct(:,:);
p=p(:,:);

zi=size(s,1);
yixi=size(s,2);
refine_ints=100;

inds=1:yixi;
fr=true(1,yixi);
delta=1e-9;

pns = nan(1,yixi);
sns = nan(1,yixi);
ctns = nan(1,yixi);

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
        ct0_stacked=ct0(ii);
        p0_stacked=p0(ii);
    end
    
    s0_stacked=s0_stacked(:,fr);
    ct0_stacked=ct0_stacked(:,fr);
    p0_stacked=p0_stacked(:,fr);
    
    pmid=0.5*(p0_stacked+p);
    bottle=gsw_rho(s0_stacked,ct0_stacked,pmid);

    cast=gsw_rho(s(:,:),ct(:,:),pmid); % 3-d density referenced to pmid
    F=cast-bottle; 
   
    [s,ct,p,sns,ctns,pns, inds,fr, dobreak]=root_core(F,delta,stack,inds,refine_ints,s,ct,p,sns,ctns,pns);
    
    if dobreak;
        break
    end
   
end


