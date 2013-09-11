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

[longs,lats] = meshgrid(longs,lats);
[e2t,e1t] = scale_fac(lats, longs);

[zi,yi,xi] = size(sa);

rho = gsw_rho(sa,ct,p_r*ones(size(sa)))-1e3;

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
    ip=find(p(:,1,1)>=1e3,1,'first');
    ptarget=p(ip,1,1)
    sns_Iak=sa(ip,:,:);
    ctns_Iak=ct(ip,:,:);
    pns_Iak=p(ip,:,:);
else
    glevels_Iak=glevels(Iak);
    % calculate properties of density surface
    if (glevels_Iak < nanmin(rho(:)))
        errordlg('Selected density surface is too light','Error Message','modal');
    elseif (glevels_Iak > nanmax(rho(:)))
        errordlg('Selected density surface is too dense','Error Message','modal');
    else
        [sns_Iak,ctns_Iak,pns_Iak,dsns_Iak,dctns_Iak,dpns_Iak] = ns_3d(sa,ct,p,rho,glevels_Iak);
    end
end

%dbstop in mld at 50
mld(1,1:yi,1:xi) = mld(sa,ct,p);

sns_Iak(pns_Iak<=mld)=nan;
ctns_Iak(pns_Iak<=mld)=nan;
pns_Iak(pns_Iak<=mld)=nan;

sns(Iak,:,:)=sns_Iak;
ctns(Iak,:,:)=ctns_Iak;
pns(Iak,:,:)=pns_Iak;

display('optimizing density surface');
tic
%dbstop in optimize_surface_exact at 624
%dbstop in optimize_surface_exact at 225
%dbstop in  optimize_surface_exact at 232 if ii==1220
%dbstop in optimize_surface_exact at 177
%dbstop in optimize_surface_exact at 92
%dbstop in  optimize_surface_exact at 152 if ii==1418
%dbstop in optimize_surface_exact at 103
sns=squeeze(sns(Iak,:,:));
ctns=squeeze(ctns(Iak,:,:));
pns=squeeze(pns(Iak,:,:));

[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:)] = optimize_surface_exact(sa,ct,p,sns,ctns,pns,e1t,e2t);
display(['optimizing density surface took ',num2str(toc),' seconds']);

        
  
