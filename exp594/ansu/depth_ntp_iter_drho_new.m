function [sns,ctns,pns] = depth_ntp_iter_drho_new(s0,ct0,p0,s,ct,p,drho)

%warning('no check of input dimensions')

refine_ints=2;

yixi=size(s,2);

inds=1:yixi;

pns = nan(1,yixi);
sns = nan(1,yixi);
ctns = nan(1,yixi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% discard land; inds contains the indices of remaining points
[s0,ct0,p0,s,ct,p,drho,inds]=discard_land(s0,ct0,p0,s,ct,p,drho,inds);

zi=size(s,1);
nxy=size(s,2);

p0_stacked=repmat(p0(:)',[zi 1]);
%keyboard
ku=sum(p0_stacked>=p,1); % adjacent bottle looking up
bottom=(ku==zi); % these are bottom bottles

i2d=1:nxy;
ku3d=ku+zi*(0:nxy-1); % 3-d

k_initial=nan*ones(1,nxy); % these are upper indices for the "bisection" (or n-section, n=refine_ints)
k3d_initial=nan*ones(1,nxy); % these are upper indices for the "bisection" (or n-section, n=refine_ints)

Fmat=nan*ones(2,nxy);

for ii=1:2*zi-1
    if ii==1
        
        % looking up
        pmid=0.5*(p0(i2d)+p(ku3d)); 
        bottle = gsw_rho(s0(i2d),ct0(i2d),pmid);
        F1=gsw_rho(s(ku3d),ct(ku3d),pmid)-bottle+drho(i2d);

        % looking down
        kd3d=ku3d+1; 
        kd3d(bottom)=1; % dummy index; remove later
        pmid=0.5*(p0(i2d)+p(kd3d)); 
        bottle = gsw_rho(s0(i2d),ct0(i2d),pmid);
        F2=gsw_rho(s(kd3d),ct(kd3d),pmid)-bottle+drho(i2d);
        F2(bottom)=nan; 
        
    else
        
        pmid=0.5*(p0(i2d)+p(ku3d)); 
        bottle = gsw_rho(s0(i2d),ct0(i2d),pmid);
        F_new=gsw_rho(s(ku3d),ct(ku3d),pmid)-bottle+drho(i2d);
        F1=F_new;
        F1(too_deep)=F_new(too_deep);
        F1(too_shallow)=F_old(too_shallow);
        F2=F_new;
        F2(too_deep)=F_old(too_deep);
        F2(too_shallow)=F_new(too_shallow);        
         
    end
    if ii==1
        moved_down=false(1,nxy);
    else
        moved_down=too_shallow;
    end
    too_deep = F1>0;
    too_shallow = F2<=0;
    inan= isnan(F1) | isnan(F2);
    todo = (too_deep | too_shallow) & ~inan;
    found_new = ~too_deep & ~too_shallow & ~inan; % F1<=0 & F2>0 & ~inan

    %keyboard 
%     if ismember(206, i2d(found_new & ~moved_down))
%         keyboard
%     end
    new_upper=found_new & ~moved_down;
    new_lower=found_new & moved_down;

    k_initial(i2d(new_upper))=ku(new_upper);
    k_initial(i2d(new_lower))=ku(new_lower)-1;
    k3d_initial(i2d(new_upper))=ku3d(new_upper);
    k3d_initial(i2d(new_lower))=ku3d(new_lower)-1;
    
    Fmat(1,i2d(found_new))=F1(found_new);
    Fmat(2,i2d(found_new))=F2(found_new);
    
    if all(~todo)
        disp('***done***')
        break
    end
    
%     if ii==2
%         break
%     end
   
    nxy_new = sum(todo);
    
    F_old=nan*ones(1,nxy_new);
    
    F_old(too_deep(todo))=F1(too_deep & todo); % F>0
    F_old(too_shallow(todo))=F2(too_shallow & todo); % F<=0  
    

    too_deep=too_deep(todo);
    too_shallow=too_shallow(todo);
    ku=ku(todo);
    ku3d=ku3d(todo);
    i2d=i2d(todo);
     
    % increment the depth index
    ku(too_deep)=ku(too_deep)-1;
    if ii==1
        ku(too_shallow)=ku(too_shallow)+2;
    else
        ku(too_shallow)=ku(too_shallow)+1;
    end
    ku3d(too_deep)=ku3d(too_deep)-1; % 3-d
    if ii==1
        ku3d(too_shallow)=ku3d(too_shallow)+2; % 3-d
    else
        ku3d(too_shallow)=ku3d(too_shallow)+1; % 3-d
    end
    
    % check if new indices are ...
    over=  ku<1; % ... above surface
    under= ku>zi; % ... below model domain
    interior= ~over & ~under;
    interior_nans=false(1,length(ku));
    interior_nans(interior)=isnan(s(ku3d(interior))); % ... in the interior but nans (bottom or sea ice)
%    keyboard
    bad = over | under | interior_nans; 

    too_deep=too_deep(~bad);
    too_shallow=too_shallow(~bad);    
    ku=ku(~bad);
    ku3d=ku3d(~bad);
    i2d=i2d(~bad);
    F_old=F_old(~bad);
    
end

%keyboard
istable=~isnan(k_initial); % stable zero crossing found
%ki=k_initial(istable); 
ki3d=k3d_initial(istable);

s0=s0(istable);
ct0=ct0(istable);
p0=p0(istable);
drho=drho(istable);
s=[s(ki3d);s(ki3d+1)];
ct=[ct(ki3d);ct(ki3d+1)];
p=[p(ki3d);p(ki3d+1)];
inds=inds(istable);
fr=true(1,sum(istable));


Fmat=Fmat(:,istable);
%disp('***end***')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


s0_stacked=repmat(s0(fr),[2 1]); % stack vertically
ct0_stacked=repmat(ct0(fr),[2 1]); 
p0_stacked=repmat(p0(fr),[2 1]);

drho_stacked=repmat(drho(fr),[2 1]);

cnt=0;
while 1
    cnt=cnt+1;
    
    if cnt~=1 % initial F values are already available
        pmid=0.5*(p0_stacked+p);
        bottle=gsw_rho(s0_stacked,ct0_stacked,pmid);

        cast=gsw_rho(s(:,:),ct(:,:),pmid); % 3-d density referenced to pmid
        F=cast-bottle+drho_stacked; 
    else
        F=Fmat;
    end
    %keyboard
    [s,ct,p,sns,ctns,pns, inds,fr]=root_core(F,inds,refine_ints,s,ct,p,sns,ctns,pns);
    
    if all(~fr) % break out of loop if all roots have been found
        break
    end
    
    s0_stacked=repmat(s0_stacked(1,fr),[refine_ints+1,1]);
    ct0_stacked=repmat(ct0_stacked(1,fr),[refine_ints+1,1]);
    p0_stacked=repmat(p0_stacked(1,fr),[refine_ints+1,1]);  
    
    drho_stacked=repmat(drho_stacked(1,fr),[refine_ints+1,1]);
    
   
end
end




function [s0,ct0,p0,s,ct,p,drho,inds]=discard_land(s0,ct0,p0,s,ct,p,drho,inds)
    iwet=~isnan(s0);

    s0=s0(iwet);
    ct0=ct0(iwet);
    p0=p0(iwet);
    drho=drho(iwet);
    s=s(:,iwet);
    ct=ct(:,iwet);
    p=p(:,iwet);
    inds=inds(iwet);
end

%         if ii==1
%             pmid=0.5*(p0_+p_);
%             bottle=gsw_rho(s0_,ct0_,pmid);
%             cast=gsw_rho(s_,ct_,pmid); % 3-d density referenced to pmid
%             F=cast-bottle+drho_; 
% 
%             F(2,bottom)=1e10; % fix bottom
% 
%             too_deep = F(1,:)>0;
%             % catch F(2,:)==0;?
%             too_shallow = F(2,:)<=0;
%             found_new = ~too_deep & ~too_shallow;
%             
%             
%         else


% s0_=repmat(s0,[2 1]); % stack vertically
% ct0_=repmat(ct0,[2 1]); 
% p0_=repmat(p0,[2 1]);
% drho_=repmat(drho,[2 1]);
% ku3d(bottom)=1; % dummy
% s_=[s(ku3d);s(ku3d+1)];
% ct_=[ct(ku3d);ct(ku3d+1)];
% p_=[p(ku3d);p(ku3d+1)];
% 
% % bottom.
% s_(1,bottom)=s(ku3d(bottom));
% ct_(1,bottom)=ct(ku3d(bottom));
% s_(2,bottom)=nan;
% ct_(2,bottom)=nan;