function [sns,ctns,pns] = depth_ntp_iter_drho_new(s0,ct0,p0,s,ct,p,drho)

%warning('no check of input dimensions')

zi=size(s,1);
yixi=size(s,2);
refine_ints=100;

inds=1:yixi;
fr=true(1,yixi);

pns = nan(1,yixi);
sns = nan(1,yixi);
ctns = nan(1,yixi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% discard land
nn=~isnan(s(:,:));
iwet=~(sum(nn,1)==0);

s0=s0(iwet);
ct0=ct0(iwet);
p0=p0(iwet);
s=s(:,iwet);
ct=ct(:,iwet);
p=p(:,iwet);
%inds=inds(iwet);
%fr=fr(iwet);


p0_stacked=repmat(p0(:)',[zi 1]);
keyboard
ku=sum(p0_stacked>=p,1); % adjacent bottle looking up
bottom=(ku==zi); % these are bottom bottles

i2d=1:sum(iwet);
ku3d=ku+zi*[0:sum(iwet)-1]; % 3-d
i2d_=i2d;
ku3d_=ku3d;

found=false(1,sum(iwet));
over=false(1,sum(iwet)); % searched up to surface
under=false(1,sum(iwet)); % searched down to bottom
stop=false(1,sum(iwet)); % no stable root

k_initial=nan*ones(1,sum(iwet));

search=true;
for ii=1:zi
    if search
        % looking up
        pmid=0.5*(p0(i2d_)+p(ku3d_)); 
        bottle = gsw_rho(s0(i2d_),ct0(i2d_),pmid);
        Fu=gsw_rho(s(ku3d_),ct(ku3d_),pmid)-bottle;

        % looking down
        kd3d=ku3d_+1; 
        kd3d(bottom)=1 % dummy index; remove later
        pmid=0.5*(p0(i2d_)+p(kd3d)); 
        Fd=gsw_rho(s(kd3d),ct(kd3d),pmid)-bottle
        Fd(bottom)=nan;

        found_new = ((Fu<=0) & (Fd>0)); % bottom bottles are lost
        k_initial(found_new)=ku(found_new);

        found(i2d_(found_new))=true;
    end
    
    %ku=ku(~found);
    inc=+ii*(-1+2*mod(ii,2));
    
    ku=ku+inc; % check for bottom or surface (or sea ice)
    over_old=over;
    under_old=under;
    error('cont here')
    over=  ku<1 | isnan(s(ku3d)); 
    under= ku>=zi;
    stop=((over & under_old) | (under & over_old)); % no stable crossing found
    
    if all(found | stop)
        break
    end
    
    eval_F = ~found & ~over & ~under & ~stop;
    if any(eval_F)
        search=true
        i2d_=i2d(eval_F)
        ku3d=ku3d+inc;
        ku3d_=ku3d(i2d_);
    else
        search=false
    end
end

disp('***end***')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% s0_stacked=repmat(s0(fr),[zi 1]); % stack vertically
% ct0_stacked=repmat(ct0(fr),[zi 1]); 
% p0_stacked=repmat(p0(fr),[zi 1]);
% 
% drho_stacked=repmat(drho(fr),[zi 1]);
% 
% cnt=0;
% while 1
%     cnt=cnt+1;
%     
%     pmid=0.5*(p0_stacked+p);
%     bottle=gsw_rho(s0_stacked,ct0_stacked,pmid);
% 
%     cast=gsw_rho(s(:,:),ct(:,:),pmid); % 3-d density referenced to pmid
%     F=cast-bottle+drho_stacked; 
%    
%     [s,ct,p,sns,ctns,pns, inds,fr]=root_core(F,inds,refine_ints,s,ct,p,sns,ctns,pns);
%     
%     if all(~fr) % break out of loop if all roots have been found
%         break
%     end
%     
%     s0_stacked=s0_stacked(1:refine_ints+1,fr);
%     ct0_stacked=ct0_stacked(1:refine_ints+1,fr);
%     p0_stacked=p0_stacked(1:refine_ints+1,fr);
%     
%     drho_stacked=drho_stacked(1:refine_ints+1,fr);
%    
% end
% 
% 
