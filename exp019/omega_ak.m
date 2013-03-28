clear all
close all
density = 'potdens';
helicity = 0;
slope = 0;
e_hel = 0;
velocity = 0;
bfrq = 0;
transport = 0;

clear ocean

load('data/input_data.mat')
lats=lats(1,:,1); longs=squeeze(longs(1,1,:))';
sa=s; clear s; % the _subs_ data saves sa in variable 's'

settings.grid = 'bgrid';
settings.wrap = 'long';
settings.slope = 'epsilon'; % density gradient error
settings.nit = 21;

glevels=[27.25]
clear g

[longs,lats] = meshgrid(longs,lats);
[e2t,e1t] = scale_fac(lats, longs,settings.wrap);

g = grav(lats);

[n2,trash] = bfrq_ct(sa,ct,p,g);

%rho = rho_from_ct(sa,ct,zeros(size(sa))) - 1e3;
rho = gpoly16ct(sa,ct); % note: sa instead of s

[zi,yi,xi] = size(sa);

%mld(1,1:yi,1:xi) = mld(sa,ct,p); % note: sa instead of s

sns = nan(length(glevels),yi,xi);
ctns = nan(length(glevels),yi,xi);
pns = nan(length(glevels),yi,xi);
dsns = nan(length(glevels),yi,xi);
dctns = nan(length(glevels),yi,xi);
dpns = nan(length(glevels),yi,xi);
sns_i = nan(length(glevels),yi,xi);
ctns_i = nan(length(glevels),yi,xi);
pns_i = nan(length(glevels),yi,xi);

Iak = 1;
% calculate properties of density surface
if (glevels(Iak) < nanmin(rho(:)))
    errordlg('Selected density surface is too light','Error Message','modal');
elseif (glevels(Iak) > nanmax(rho(:)))
    errordlg('Selected density surface is too dense','Error Message','modal');
else
    [sns(Iak,:,:),ctns(Iak,:,:),pns(Iak,:,:),dsns(Iak,:,:),dctns(Iak,:,:),dpns(Iak,:,:)] = ns_3d(sa,ct,p,rho,glevels(Iak));
end

mld(1,1:yi,1:xi) = mld(sa,ct,p);
[sns(Iak,:,:)] = cut_off(squeeze(sns(Iak,:,:)),squeeze(pns(Iak,:,:)),squeeze(mld));
[ctns(Iak,:,:)] = cut_off(squeeze(ctns(Iak,:,:)),squeeze(pns(Iak,:,:)),squeeze(mld));
[pns(Iak,:,:)] = cut_off(squeeze(pns(Iak,:,:)),squeeze(pns(Iak,:,:)),squeeze(mld));

display('optimizing density surface');
tic
%[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:),ithist] = optimize_surface_exact(sa,ct,p,g,n2,sns(Iak,:,:),ctns(Iak,:,:),pns(Iak,:,:),e1t,e2t,settings);
[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:),ithist] = optimize_surface_exact(sa,ct,p,g,n2,sns(Iak,:,:),ctns(Iak,:,:),pns(Iak,:,:),e1t,e2t,settings);
display(['optimizing density surface took ',int2str(toc),' seconds']);

% save variables for postprocessing 
save('data/ansu_hist.mat', 'ithist')
        
  
