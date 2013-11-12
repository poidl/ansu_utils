function [sns_out,ctns_out,pns_out] = dz_from_drho(sns, ctns, pns, s, ct, p, drho );

[zi,yi,xi]=size(s);
drho = permute(drho, [3 1 2]);

delta = 1e-11;

rho_surf=gsw_rho(sns(:),ctns(:),pns(:));
t2=rho_surf-drho(:);

inds=1:yi*xi;
fr=true(1,yi*xi);
pns_out = nan(yi,xi);
sns_out = nan(yi,xi);
ctns_out = nan(yi,xi);
pns=pns(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative code

% stack=zi;
% ii=bsxfun(@times,1:yi*xi,ones(stack,1));
% pns_stacked=pns(ii); % stack pressure of current surface vertically
% t2_stacked_full=t2(ii); % stack locally referenced density of current surface vertically
% 
% t1=gsw_rho(s(:,:),ct(:,:),pns_stacked(:,inds)); % 3-d density referenced to pressure of the current surface
% t2_3d=t2_stacked_full(:,inds); %
% tni=t1-t2_3d; % rho-(rho_s+rho'); find corrected surface by finding roots of this term
% 
% Itni_p = tni>0;
% Itni_n = tni<0;
% 
% zc = any(Itni_p,1) & any(Itni_n,1); % horizontal indices of locations where zero-crossing occurs
% [min_tni, lminr]=min(abs(tni));
% cond1=min_tni>delta;
% final=min_tni<=delta; % at these horizontal indices root has been found
% fr= zc & cond1; %  at these horizontal locations we have to use root finding
% 
% lminr=lminr+stack*[0:size(Itni_n,2)-1];
% lminr=lminr(final);
% 
% sns_out(inds(final))=s(lminr); % adjust surface where root has already been found
% ctns_out(inds(final)) =ct(lminr);
% pns_out(inds(final)) =p(lminr);
% 
% k=sum(Itni_n,1); % find indices of flattened 3d-array where we use root finding
% k=k+stack*[0:size(Itni_n,2)-1];
%
% options=optimset('TolX',1e-3);
% for ii=inds(fr);
%     su=s(k(ii)); sl=s(k(ii)+1); % linear interpolation of s and ct
%     ctu=ct(k(ii)); ctl=ct(k(ii)+1);
%     pu=p(k(ii)); pl=p(k(ii)+1);
%     
%     sp=@(p) su+(sl-su)*(p-pu)/(pl-pu);
%     ctp=@(p) ctu+(ctl-ctu)*(p-pu)/(pl-pu);
%     fac=1./sqrt(gsw_rho(su,ctu,pu)^2+gsw_rho(sl,ctl,pl)^2); % scale to avoid having to set tolerance ?
%     drho_normalized=@(p)  fac.*(gsw_rho(sp(p), ctp(p), pns(ii)) -t2(ii));
%     
%     %dmyrho_normalized=@(p,su,sl,ctu,ctl,pu,pl)  1.*(gsw_rho( su+(sl-su)*(p-pu)/(pl-pu), ctu+(ctl-ctu)*(p-pu)/(pl-pu), pns(ii))-t2(ii));
%     %drho_normalized=@(p) dmyrho_normalized(p,su,sl,ctu,ctl,pu,pl);
%     proot=fzero(drho_normalized, [pu,pl],options);
%     
%     sns_out(ii)= su+(sl-su)*(proot-pu)/(pl-pu);
%     ctns_out(ii)=ctu+(ctl-ctu)*(proot-pu)/(pl-pu);
%     pns_out(ii)=proot;
% end

% alternative code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative code
    
refine_ints=100;

cnt=0;
while 1
    cnt=cnt+1;
    
    if cnt==1 % in first iteration pns_l is stacked vertically zi times, after that it is stacked refine_ints times
        stack=zi;
    elseif cnt==2
        stack=refine_ints+1;
    end
    if cnt==1 | cnt==2
        ii=bsxfun(@times,1:yi*xi,ones(stack,1));
        pns_stacked=pns(ii); % stack pressure of current surface vertically
        t2_stacked=t2(ii); % stack locally referenced density of current surface vertically
    end
    
    pns_stacked=pns_stacked(:,fr);
    t2_stacked=t2_stacked(:,fr);
    
    t1=gsw_rho(s(:,:),ct(:,:),pns_stacked); % 3-d density referenced to pressure of the current surface

    F=t1-t2_stacked; % rho-(rho_s+rho'); find corrected surface by finding roots of this term
    
    %dbstop in root_core at 11
    [s,ct,p,sns_out,ctns_out,pns_out, inds, fr, dobreak]=root_core(F,delta,stack,inds,refine_ints,s,ct,p,sns_out,ctns_out,pns_out);
    
    if dobreak;
        break
    end

end

% alternative code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
