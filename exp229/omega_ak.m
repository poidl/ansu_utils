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

% [zi,yi,xi]=size(s);
% 
% ii=yi*45+10;
% ii=[ii:ii+2];
% inds=[ ii ii+yi ii+2*yi ];
% s=s(:,inds);
% ct=ct(:,inds);
% p=p(:,inds);
% s=reshape(s,[size(s,1), 3 3]);
% ct=reshape(ct,[size(s,1), 3 3]);
% p=reshape(p,[size(s,1), 3 3]);

% s(:,5)=nan;
% ct(:,5)=nan;
% p(:,5)=nan;
 %s(:,:,2)=nan;
 %ct(:,:,2)=nan;
 %p(:,:,2)=nan;
% s(:,2,:)=nan;
% ct(:,2,:)=nan;
% p(:,2,:)=nan;


% keep=[[yi*(xi-1)+1:yi*xi]'  [1:yi]'];
% setnan=true(1,xi*yi);
% setnan(keep)=false;
% s(:,setnan)=nan;
% ct(:,setnan)=nan;
% p(:,setnan)=nan;
% ct=ct(:,15:20,:);
% p=p(:,15:20,:);
% lats=lats(:,15:20,:);
% longs=longs(:,15:20,:);

% mydims=[size(s,1), size(keep)];
% s=reshape(s,mydims);
% ct=reshape(ct,mydims);
% p=reshape(p,mydims);
% lats=lats(:,keep);
% longs=longs(:,keep);
% reg=regions{1};
% reg= intersect(reg, keep);
% regions{1}=reg;

lats=lats(1,:,1); longs=squeeze(longs(1,1,:))';
sa=s; clear s; % the _subs_ data saves sa in variable 's'

[longs,lats] = meshgrid(longs,lats);

[zi,yi,xi] = size(sa);


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
    rho = gsw_rho(sa,ct,p_r*ones(size(sa)))-1e3;
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

mld(1,1:yi,1:xi) = mld(sa,ct,p);

sns_Iak(pns_Iak<=mld)=nan;
ctns_Iak(pns_Iak<=mld)=nan;
pns_Iak(pns_Iak<=mld)=nan;

sns(Iak,:,:)=sns_Iak;
ctns(Iak,:,:)=ctns_Iak;
pns(Iak,:,:)=pns_Iak;

display('optimizing density surface');
tic

%dbstop in  optimize_surface_exact at 232 if ii==1220
%dbstop in optimize_surface_exact at 108
dbstop in optimize_surface_exact at 125

sns=squeeze(sns(Iak,:,:));
ctns=squeeze(ctns(Iak,:,:));
pns=squeeze(pns(Iak,:,:));

save_netcdf(sa,ct,p,sns,ctns,pns);
[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:)] = optimize_surface_exact(sa,ct,p,sns,ctns,pns);
display(['optimizing density surface took ',num2str(toc),' seconds']);

        
  
