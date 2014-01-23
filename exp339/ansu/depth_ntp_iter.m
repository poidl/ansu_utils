function [sns,ctns,pns] = depth_ntp_iter(s0,ct0,p0,s,ct,p)

%warning('no check of input dimensions')

zi=size(s,1);
yixi=size(s,2);
refine_ints=100;

inds=1:yixi;
fr=true(1,yixi);

pns = nan(1,yixi);
sns = nan(1,yixi);
ctns = nan(1,yixi);

s0_stacked=repmat(s0(fr),[zi 1]); % stack vertically
ct0_stacked=repmat(ct0(fr),[zi 1]); 
p0_stacked=repmat(p0(fr),[zi 1]);

cnt=0;
while 1
    cnt=cnt+1;
    
    pmid=0.5*(p0_stacked+p);
    bottle=gsw_rho(s0_stacked,ct0_stacked,pmid);

    cast=gsw_rho(s(:,:),ct(:,:),pmid); % 3-d density referenced to pmid
    F=cast-bottle; 
   
    [s,ct,p,sns,ctns,pns, inds,fr]=root_core(F,inds,refine_ints,s,ct,p,sns,ctns,pns);
    
    if all(~fr) % break out of loop if all roots have been found
        break
    end

    s0=s0(fr);
    ct0=ct0(fr);
    p0=p0(fr);
    drho=drho(fr);
    
    s0_stacked=repmat(s0,[refine_ints+1 1]); % stack vertically
    ct0_stacked=repmat(ct0,[refine_ints+1 1]); 
    p0_stacked=repmat(p0,[refine_ints+1 1]);

   
end


