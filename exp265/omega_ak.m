clear all
close all

params; % load some parameters

load('data/input_data.mat')
for ii=1:size(s,3)
    for jj=1:size(s,2)
        kk=find(isnan(s(:,jj,ii)),1,'first');
        s(kk:end,jj,ii)=nan;
        ct(kk:end,jj,ii)=nan;
%        p(kk:end,jj,ii)=nan;
    end
end

lats=lats(1,:,1); longs=squeeze(longs(1,1,:))';
sa=s; clear s; % the _subs_ data saves sa in variable 's'

load('../exp258/data/iteration_history.mat','sns_hist');
regions=find_regions(squeeze(sns_hist(end,:,:)));

[zi,yi,xi]=size(sa);
imaxregion=0; % index of largest region
npmaxreg=0; % number of points in largest region
for i=1:length(regions)
    if length(regions{i})>npmaxreg;
        npmaxreg=length(regions{i});
        imaxreg=i;
    end
end
setnan=true(1,xi*yi);
setnan(regions{imaxreg})=false;
s(:,setnan)=nan;
ct(:,setnan)=nan;
p(:,setnan)=nan;
               
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % helicity spike
% ilon=find(longs>=180,1,'first');
% ilat=find(lats>=-60,1,'first');
% 
% ct_cast=squeeze(ct(:,ilat,ilon));
% s_cast=squeeze(sa(:,ilat,ilon));
% r=gsw_rho_CT_exact(s_cast, ct_cast,p(:,ilat,ilon));
% 
% ctd=linspace(10.0, 2.0, size(ct,1));
% ctnew_cast=ct_cast+ctd';
% snew_cast=gsw_SA_from_rho_CT_exact(r,ctnew_cast,p(:,ilat,ilon));
% 
% % check
% rnew=gsw_rho_CT_exact(snew_cast, ctnew_cast,p(:,ilat,ilon));
% disp(['Max. dens. deviation: ', num2str(max(abs(r-rnew)))])
% 
% ct(:,ilat,ilon)=ctnew_cast;
% sa(:,ilat,ilon)=snew_cast;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sns = nan(length(nlevels),yi,xi);
ctns = nan(length(nlevels),yi,xi);
pns = nan(length(nlevels),yi,xi);
dsns = nan(length(nlevels),yi,xi);
dctns = nan(length(nlevels),yi,xi);
dpns = nan(length(nlevels),yi,xi);
sns_i = nan(length(nlevels),yi,xi);
ctns_i = nan(length(nlevels),yi,xi);
pns_i = nan(length(nlevels),yi,xi);

Iak = 1;

if initial_surface_at_constant_pressure;
    ip=find(p(:,1,1)>=initial_pressure,1,'first');
    ptarget=p(ip,1,1)
    sns_Iak=sa(ip,:,:);
    ctns_Iak=ct(ip,:,:);
    pns_Iak=p(ip,:,:);
    pns_Iak(isnan(sns_Iak))=nan;
else
    glevels_Iak=1e3+glevels(Iak);
    delta=1e-11;
    [zi,yi,xi]=size(sa);
    
    s_2=sa(:,:);
    ct_2=ct(:,:);
    p_2=p(:,:);
    inds=1:yi*xi;
    fr=true(1,yi*xi);
    pns_out = nan(yi,xi);
    sns_out = nan(yi,xi);
    ctns_out = nan(yi,xi);
    
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
        end
        
        F=gsw_rho(s_2,ct_2,p_r*ones(size(s_2)))-glevels_Iak;
        [s_2,ct_2,p_2,sns_out,ctns_out,pns_out, inds, fr, dobreak]=root_core(F,delta,stack,inds,refine_ints,s_2,ct_2,p_2,sns_out,ctns_out,pns_out);
        
        if dobreak;
            break
        end
    end
end

sns=sns_out;
ctns=ctns_out;
pns=pns_out;

display('optimizing density surface');
tic
%dbstop in optimize_surface_exact at 624
%dbstop in optimize_surface_exact at 225
%dbstop in  optimize_surface_exact at 232 if ii==1220
%dbstop in optimize_surface_exact at 177
%dbstop in optimize_surface_exact at 104
%dbstop in optimize_surface_exact at 443

load('../exp259/data/iteration_history.mat','sns_hist','pns_hist','ctns_hist');
sns=squeeze(sns_hist(end,:,:));
ctns=squeeze(ctns_hist(end,:,:));
pns=squeeze(pns_hist(end,:,:));
sns(setnan)=nan;
ctns(setnan)=nan;
pns(setnan)=nan;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% localized perturbation
drho=zeros(size(sns));
drho(10,:)=0.2;
drho(11,:)=-0.2;
drho(30,:)=0.2;
drho(31,:)=-0.2;
[sns,ctns,pns] = dz_from_drho(sns, ctns, pns, sa, ct, p, drho );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:)] = optimize_surface_exact(sa,ct,p,sns,ctns,pns);
display(['optimizing density surface took ',num2str(toc),' seconds']);

        
        
  
