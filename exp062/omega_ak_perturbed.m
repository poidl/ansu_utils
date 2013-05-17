clear all
close all

user_input; % load configuration from settings.m

load('data/input_data_perturbed.mat')
load('data/surface_optimized')
lats=lats(1,:,1); longs=squeeze(longs(1,1,:))';
sa=s; clear s; % the _subs_ data saves sa in variable 's'

[longs,lats] = meshgrid(longs,lats);
[e2t,e1t] = scale_fac(lats, longs);

g = gsw_grav(lats);

[zi,yi,xi] = size(sa);

lats=repmat(permute(lats,[3,1,2]),[zi,1,1]);
[n2tmp,p_mid] = gsw_Nsquared(sa,ct,p,lats);
n2=reshape(n2tmp,size(sa)-[1 0 0]);

Iak = 1;

sns(Iak,:,:)=sns_i;
ctns(Iak,:,:)=ctns_i;
pns(Iak,:,:)=pns_i;


sns_i = nan(length(glevels),yi,xi);
ctns_i = nan(length(glevels),yi,xi);
pns_i = nan(length(glevels),yi,xi);

display('optimizing density surface');
tic
%dbstop in optimize_surface_exact at 624
%dbstop in optimize_surface_exact at 69
[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:)] = optimize_surface_exact(sa,ct,p,g,n2,sns(Iak,:,:),ctns(Iak,:,:),pns(Iak,:,:),e1t,e2t);
display(['optimizing density surface took ',num2str(toc),' seconds']);

        
  
