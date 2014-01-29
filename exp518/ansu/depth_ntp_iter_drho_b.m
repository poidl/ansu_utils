function [sns,ctns,pns] = depth_ntp_iter_drho_b(s0,ct0,p0,s,ct,p,drho,b,pb)

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

drho_stacked=repmat(drho(fr),[zi 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myb: interpolate b on p
    myb=nan*ones(size(s));
    for kk=1:size(s,1)
        p_slice=p(kk,:); p_slice=repmat(p_slice, [size(b,1) 1]);
        up=pb<p_slice;
        kup=sum(up,1);
       
        inan=isnan(p(kk,:));
        tooshallow=(kup==0 & ~inan);
        toodeep=(kup==size(pb,1) & ~inan);
        bad=(tooshallow | toodeep | inan); % dummy index; remove later
        kup(bad)=1; % dummy index; remove later
        
        kup=kup+size(b,1)*[0:size(b,2)-1]; % 3d
        
        b_1=b(kup);
        b_2=b(kup+1);
        pb_1=pb(kup);
        pb_2=pb(kup+1);
        
        dp=(p(kk,:)-pb_1)./(pb_2-pb_1);
        dp(bad)=nan;
        myb(kk,:)=b_1+(b_2-b_1).*dp;
        
        myb(kk,inan)=nan;
        myb(kk,tooshallow)=b(1,tooshallow);
        myb(kk,toodeep)=b(end,toodeep);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% sb: interpolate b onto surface

    p_slice=repmat(p0, [size(b,1) 1]);
    up=pb<p_slice;
    kup=sum(up,1);

    inan=isnan(p0);
    tooshallow=(kup==0 & ~inan);
    bad=(tooshallow | inan); % dummy index; remove later
    kup(bad)=1; % dummy index; remove later

    kup=kup+size(b,1)*[0:size(b,2)-1]; % 3d

    b_1=b(kup);
    b_2=b(kup+1);
    pb_1=pb(kup);
    pb_2=pb(kup+1);

    dp=(p0-pb_1)./(pb_2-pb_1);
    dp(bad)=nan;
    sb=b_1+(b_2-b_1).*dp;

    sb(inan)=nan;
    sb(tooshallow)=b(1,tooshallow); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    sb_stacked=repmat(sb(fr), [zi 1]); 

cnt=0;
while 1
    cnt=cnt+1;
    
    pmid=0.5*(p0_stacked+p);
    bottle=gsw_rho(s0_stacked,ct0_stacked,pmid);

    cast=gsw_rho(s(:,:),ct(:,:),pmid); % 3-d density referenced to pmid
        
    F=0.5*(sb_stacked+myb).*(cast-bottle)+drho_stacked;  
    %F=sb_stacked.*(cast-bottle)+drho_stacked;
    %F=cast-bottle+drho_stacked./sb_stacked;
    
    [s,ct,p,myb,sns,ctns,pns, inds,fr]=root_core_b(F,inds,refine_ints,s,ct,p,myb,sns,ctns,pns);
    
    if all(~fr) % break out of loop if all roots have been found
        break
    end
    
    s0=s0(fr);
    ct0=ct0(fr);
    p0=p0(fr);
    drho=drho(fr);
    sb=sb(fr);
    
    s0_stacked=repmat(s0,[refine_ints+1 1]); % stack vertically
    ct0_stacked=repmat(ct0,[refine_ints+1 1]); 
    p0_stacked=repmat(p0,[refine_ints+1 1]);

    drho_stacked=repmat(drho,[refine_ints+1 1]);
    sb_stacked=repmat(sb,[refine_ints+1 1]); 
   
end


