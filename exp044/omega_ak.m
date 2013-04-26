clear all
close all
density = 'potdens';
helicity = 0;
slope = 0;
e_hel = 0;
velocity = 0;
bfrq = 0;
transport = 0;

load('data/input_data.mat')
lats=lats(1,:,1); longs=squeeze(longs(1,1,:))';
sa=s; clear s; % the _subs_ data saves sa in variable 's'

settings.grid = 'bgrid';
settings.wrap = 'long';
settings.slope = 'epsilon'; % density gradient error
settings.solver = 'iterative';
settings.nit = 21;

glevels=[27.5]

[longs,lats] = meshgrid(longs,lats);
[e2t,e1t] = scale_fac(lats, longs,settings.wrap);

g = gsw_grav(lats);

[zi,yi,xi] = size(sa);

lats=repmat(permute(lats,[3,1,2]),[zi,1,1]);
[n2tmp,p_mid] = gsw_Nsquared(sa,ct,p,lats);
n2=reshape(n2tmp,size(sa)-[1 0 0]);

rho = gsw_rho(sa,ct,zeros(size(sa)))-1e3;

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
glevels_Iak=glevels(Iak);

% calculate properties of density surface
if (glevels_Iak < nanmin(rho(:)))
    errordlg('Selected density surface is too light','Error Message','modal');
elseif (glevels_Iak > nanmax(rho(:)))
    errordlg('Selected density surface is too dense','Error Message','modal');
else
    [sns_Iak,ctns_Iak,pns_Iak,dsns_Iak,dctns_Iak,dpns_Iak] = ns_3d(sa,ct,p,rho,glevels_Iak);
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
%dbstop in optimize_surface_exact at 155
[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:),ithist] = optimize_surface_exact(sa,ct,p,g,n2,sns(Iak,:,:),ctns(Iak,:,:),pns(Iak,:,:),e1t,e2t,settings);
display(['optimizing density surface took ',num2str(toc),' seconds']);

% save variables for postprocessing 
save('data/ansu_hist.mat', 'ithist')
        
  
