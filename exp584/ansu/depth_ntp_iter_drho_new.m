function [sns,ctns,pns] = depth_ntp_iter_drho_new(s0,ct0,p0,s,ct,p,drho)

%warning('no check of input dimensions')

zi=size(s,1);
yixi=size(s,2);
refine_ints=2;

inds=1:yixi;
fr=true(1,yixi);

pns = nan(1,yixi);
sns = nan(1,yixi);
ctns = nan(1,yixi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% discard land
nn=~isnan(s0);
iwet=~(sum(nn,1)==0);

s0=s0(iwet);
ct0=ct0(iwet);
p0=p0(iwet);
drho=drho(iwet);
s=s(:,iwet);
ct=ct(:,iwet);
p=p(:,iwet);
inds=inds(iwet);
fr=fr(iwet);

p0_stacked=repmat(p0(:)',[zi 1]);
%keyboard
ku=sum(p0_stacked>=p,1); % adjacent bottle looking up
bottom=(ku==zi); % these are bottom bottles

i2d=1:sum(iwet);
ku3d=ku+zi*[0:sum(iwet)-1]; % 3-d
i2d_=i2d;
ku3d_=ku3d;

found=false(1,sum(iwet));
bad=false(1,sum(iwet));

k_initial=nan*ones(1,sum(iwet)); % these are upper indices for the "bisection" (or n-section, n=refine_ints)
k3d_initial=nan*ones(1,sum(iwet)); % these are upper indices for the "bisection" (or n-section, n=refine_ints)

Fmat=nan*ones(2,sum(iwet));

%%
%drho=0*drho;
%drho(23)=-0.165;
%%
%keyboard
%keyboard
% find k_initial
search=true;
for ii=1:zi
    if search
        % looking up
        pmid=0.5*(p0(i2d_)+p(ku3d_)); 
        bottle = gsw_rho(s0(i2d_),ct0(i2d_),pmid);
        Fu=gsw_rho(s(ku3d_),ct(ku3d_),pmid)-bottle+drho(i2d_);

        % looking down
        kd3d=ku3d_+1; 
        kd3d(bottom)=1; % dummy index; remove later
        pmid=0.5*(p0(i2d_)+p(kd3d)); 
        bottle = gsw_rho(s0(i2d_),ct0(i2d_),pmid);
        Fd=gsw_rho(s(kd3d),ct(kd3d),pmid)-bottle+drho(i2d_);
        Fd(bottom)=nan;

        found_new = ((Fu<=0) & (Fd>0)); % bottom bottles are lost
        k_initial(i2d_(found_new))=ku(i2d_(found_new));
        k3d_initial(i2d_(found_new))=ku3d(i2d_(found_new));
        found(i2d_(found_new))=true;
        Fmat(1,i2d_(found_new))=Fu(found_new);
        Fmat(2,i2d_(found_new))=Fd(found_new);
    end
    
    %ku=ku(~found);
    inc=+ii*(-1+2*mod(ii,2));
    
    % increment the depth index and check if new indices ...
    ku=ku+inc; % increment
    ku3d=ku3d+inc; % 3-d
    bad_old=bad;
    over=  ku<1; % ... are above surface
    under= ku>zi-1; % ... are below model domain
    interior= ~over & ~under;
    interior_nans=false(1,length(ku));
    interior_nans(interior)=isnan(s(ku3d(interior))); % are in the interior but nans (bottom or sea ice) TODO: check lower bottle too?
    bad = over | under | interior_nans; 
    stop= bad & bad_old; % undefined values in both directions.
    
    if all(found | stop) % either a bottle pair has been found or searching in both directions was unsuccessfully completed
        break
    end
    
    eval_F = ~found & ~bad;
    if any(eval_F)
        search=true;
        i2d_=i2d(eval_F);
        ku3d_=ku3d(eval_F);
    else
        search=false;
    end
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
fr=fr(istable);

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
    
    pmid=0.5*(p0_stacked+p);
    bottle=gsw_rho(s0_stacked,ct0_stacked,pmid);

    cast=gsw_rho(s(:,:),ct(:,:),pmid); % 3-d density referenced to pmid
    F=cast-bottle+drho_stacked; 
   
    [s,ct,p,sns,ctns,pns, inds,fr]=root_core(F,inds,refine_ints,s,ct,p,sns,ctns,pns);
    
    if all(~fr) % break out of loop if all roots have been found
        break
    end
    
    s0_stacked=repmat(s0_stacked(1,fr),[refine_ints+1,1]);
    ct0_stacked=repmat(ct0_stacked(1,fr),[refine_ints+1,1]);
    p0_stacked=repmat(p0_stacked(1,fr),[refine_ints+1,1]);  
    
    drho_stacked=repmat(drho_stacked(1,fr),[refine_ints+1,1]);
    
   
end


