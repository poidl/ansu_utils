clear all
close all
density = 'potdens';
helicity = 0;
slope = 0;
e_hel = 0;
velocity = 0;
bfrq = 0;
transport = 0;
% e_therm_cab = 0;
%nit = 180;
which_surface = 'initial';
streamfunc = 0;
% geo_vel = 0;
global ocean_options settings%#ok

clear ocean

load('data/input_data.mat')

sa=s; clear s; % the _subs_ data saves sa in variable 's'
%todo: the following has to work with all kinds of resolution!
if min(longs(:)) <= 1 & max(longs(:)) >=356
    wrap = 'long';
else
    wrap = 'none';
end

grid = 'bgrid';


settings.grid = grid;
settings.wrap = wrap;
settings.keyword = 'op';
settings.solver = 'iterative';
settings.density = 'poly'; % gamma^rf
settings.slope = 'epsilon'; % density gradient error
settings.nit = 21;
ref_level = 1800;

glevels=[27.25]
clear g

[e2t,e1t] = scale_fac(squeeze(lats(1,:,:)),squeeze(longs(1,:,:)));

g = gsw_grav(squeeze(lats(1,:,:)));

[n2,p_mid] = gsw_Nsquared(sa,ct,p,lats);
n2=reshape(n2,size(sa)-[1 0 0]);
n2_smth = nan*ones(size(n2));
for ii= 1:size(n2,2)
    for jj = 1:size(n2,3)
        n2_smth(:,ii,jj) = smooth(n2(:,ii,jj),7);
    end
end

rho = gsw_rho(sa,ct,zeros(size(sa)))-1e3;

[zi,yi,xi] = size(sa);
%todo: mld takes sp, not sa
mld(1,1:yi,1:xi) = mld(sa,ct,p);

sns = nan(length(glevels),yi,xi);
ctns = nan(length(glevels),yi,xi);
pns = nan(length(glevels),yi,xi);
dsns = nan(length(glevels),yi,xi);
dctns = nan(length(glevels),yi,xi);
dpns = nan(length(glevels),yi,xi);
sns_i = nan(length(glevels),yi,xi);
ctns_i = nan(length(glevels),yi,xi);
pns_i = nan(length(glevels),yi,xi);
n2_ns_i = nan(length(glevels),yi,xi);
e_cab_i = nan(length(glevels),yi,xi);
e_therm_i = nan(length(glevels),yi,xi);
e_cab_trans_i = nan(length(glevels),yi,xi);
e_cab_trans_sum_i = nan(length(glevels),yi,xi);
streamfunc_i = nan(length(glevels),yi,xi);
geo_vel_x_i = nan(length(glevels),yi,xi);
geo_vel_y_i = nan(length(glevels),yi,xi);

%for Iak=1:length(glevels)

%    display(['Iak',int2str(Iak)])
 %   disp(glevels(Iak))
        Iak = 1;
        % calculate properties of density surface
        if (glevels(Iak) < nanmin(rho(:)))
            errordlg('Selected density surface is too light','Error Message','modal');
        elseif (glevels(Iak) > nanmax(rho(:)))
            errordlg('Selected density surface is too dense','Error Message','modal');
        else
            [sns(Iak,:,:),ctns(Iak,:,:),pns(Iak,:,:),dsns(Iak,:,:),dctns(Iak,:,:),dpns(Iak,:,:)] = ns_3d(sa,ct,p,rho,glevels(Iak));
        end
        
        % calculate average pressure of density surface
        %p_ave = nanmean(pns(:));
        %     display(['the average pressure of this surface is ',int2str(p_ave),' meters']);
        
        [sns(Iak,:,:)] = cut_off(squeeze(sns(Iak,:,:)),squeeze(pns(Iak,:,:)),squeeze(mld));
        [ctns(Iak,:,:)] = cut_off(squeeze(ctns(Iak,:,:)),squeeze(pns(Iak,:,:)),squeeze(mld));
        [pns(Iak,:,:)] = cut_off(squeeze(pns(Iak,:,:)),squeeze(pns(Iak,:,:)),squeeze(mld));
        
        display('optimizing density surface');
        tic
        [sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:),ithist] = optimize_surface_exact(sa,ct,p,g,n2,sns(Iak,:,:),ctns(Iak,:,:),pns(Iak,:,:),e1t,e2t,settings);
        display(['optimizing density surface took ',int2str(toc),' seconds']);
        
        % save variables for postprocessing 
        save('data/ansu_hist.mat', 'ithist')
        
        if 0
            %     tic
            %     display('calculating buoyancy frequency on omega surface');
            [n2_ns_i(Iak,:,:)] = var_on_surf(pns_i(Iak,:,:),p_mid,n2);
            %     display(['calculating buoyancy frequency on omega surface took ',int2str(toc),' seconds']);

            %     tic
            %     display('calculating diapycnal velocity due to thermobaricity/cabbeling through omega surface');
            [e_cab_i(Iak,:,:),e_therm_i(Iak,:,:)] = e_therm_cab(sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:), ...
                n2_ns_i(Iak,:,:),g,e1t,e2t);
            %     display(['calculating diapycnal velocity due to thermobaricity/cabbeling through surface took ',int2str(toc),' seconds']);


            %     tic
            %     display('calculating diapycnal transport on omega surface');
            [e_cab_trans_i(Iak,:,:),e_cab_trans_sum_i(Iak,:,:)] = dia_trans(e_cab_i(Iak,:,:),e1t,e2t);
            %     display(['calculating diapycnal transport on omega surface took ',int2str(toc),' seconds']);

            %     tic
            %     display('calculating geostrophic streamfunction on omega surface');
            [streamfunc_i(Iak,:,:)] = optimize_streamfunc(s,ct,p,sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:),ref_level);
            %     display(['calculating geostrophic streamfunction on omega surface took ',int2str(toc),' seconds']);

            %     tic
            %     display('calculating geostrophic velocities on omega surface');
            sns_i_dummy = sns_i(Iak,:,:);
            ctns_i_dummy = ctns_i(Iak,:,:);
            [geo_vel_x_i(Iak,:,:),geo_vel_y_i(Iak,:,:)] = geo_vel(streamfunc_i(Iak,:,:),s,ct,p,nanmean(sns_i_dummy(:)),nanmean(ctns_i_dummy(:)),lats,e1t,e2t,ref_level);
            %     display(['calculating geostrophic velocities on omega surface took ',int2str(toc),' seconds']);

            %     tic
            %     display('calculating buoyancy frequency on  omega surface');
            %[n2_ns] = var_on_surf(pns,p_mid,n2);
            [n2_ns_i(Iak,:,:)] = var_on_surf(pns_i(Iak,:,:),p_mid,n2);
            %     display(['calculating buoyancy frequency on initial and omega surface took ',int2str(toc),' seconds']);

            %     tic
            %     display('calculating diapycnal transports through omega surface');
            [e_cab_trans_i(Iak,:,:),e_cab_trans_sum_i(Iak,:,:)] = dia_trans(e_cab_i(Iak,:,:),e1t,e2t);
            %     display(['calculating diapycnal transports through omega surface took ',int2str(toc),' seconds']);

            %save (['/home/c000573-hf/bar747/all/omega_maps/omega_',num2str(Iak),'.mat'])
        end
%end
%save (['omega_gk2.mat'])