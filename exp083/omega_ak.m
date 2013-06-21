clear all
close all

user_input; % load configuration from settings.m

load('data/input_data.mat')
lats=lats(1,:,1); longs=squeeze(longs(1,1,:))';
sa=s; clear s; % the _subs_ data saves sa in variable 's'

[longs,lats] = meshgrid(longs,lats);
[e2t,e1t] = scale_fac(lats, longs);

g = gsw_grav(lats);

[zi,yi,xi] = size(sa);

lats=repmat(permute(lats,[3,1,2]),[zi,1,1]);
[n2tmp,p_mid] = gsw_Nsquared(sa,ct,p,lats);
n2=reshape(n2tmp,size(sa)-[1 0 0]);

rho = gsw_rho(sa,ct,zeros(size(sa)))-1e3;

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
%dbstop in optimize_surface_exact at 201
[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:)] = optimize_surface_exact(sa,ct,p,g,n2,sns(Iak,:,:),ctns(Iak,:,:),pns(Iak,:,:),e1t,e2t);
display(['optimizing density surface took ',num2str(toc),' seconds']);
vars= {'sns_i','ctns_i','pns_i'};
save('data/surface_optimized.mat',vars{:})

        
  
