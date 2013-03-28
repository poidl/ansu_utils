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
global ocean_options longs lats settings%#ok

clear ocean

grid = 'bgrid';

load gk_ak_gamma.mat
if min(min(longs)) <= 1 & max(max(longs)) >=359
    wrap = 'long';
else
    wrap = 'none';
end

settings.grid = grid;
settings.wrap = wrap;
settings.keyword = 'op';
settings.solver = 'iterative';
settings.density = 'poly'; % gamma^rf
settings.slope = 'epsilon'; % density gradient error
%settings.nit = 180;
settings.nit=21;
ref_level = 1800;

longss = longs;
latss = lats;
clear longs lats
[longs,lats] = meshgrid(longss,latss);
clear latss longss
%glevels = [26:0.1:28];
glevels=[27.25]
clear g

[e2t,e1t] = scale_fac(lats,longs);

g = grav(lats);

[n2,p_mid] = bfrq_ct(s,ct,p,g);
n2_smth = nan(400,171,360);
if 1 %skip smoothing to debug quickly
    for ii= 1:171
        for jj = 1:360
            n2_smth(:,ii,jj) = smooth(n2(:,ii,jj),7);
        end
    end
else
    n2_smth=n2(1:end-1,:,:); 
end

rho = gpoly16ct(s,ct);
poly = rho;

[zi,yi,xi] = size(s);
mld(1,1:yi,1:xi) = mld(s,ct,p);

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

%for Iak = round(length(glevels)/2):round(3*(length(glevels)/4))
for Iak=1:1
    %try
        % calculate properties of density surface
        if (glevels(Iak) < nanmin(rho(:)))
            errordlg('Selected density surface is too light','Error Message','modal');
        elseif (glevels(Iak) > nanmax(rho(:)))
            errordlg('Selected density surface is too dense','Error Message','modal');
        else
            [sns(Iak,:,:),ctns(Iak,:,:),pns(Iak,:,:),dsns(Iak,:,:),dctns(Iak,:,:),dpns(Iak,:,:)] = ns_3d(s,ct,p,rho,glevels(Iak));
        end
        
        % calculate average pressure of density surface
        %p_ave = nanmean(pns(:));
        %     display(['the average pressure of this surface is ',int2str(p_ave),' meters']);
        
        [sns(Iak,:,:)] = cut_off(squeeze(sns(Iak,:,:)),squeeze(pns(Iak,:,:)),squeeze(mld));
        [ctns(Iak,:,:)] = cut_off(squeeze(ctns(Iak,:,:)),squeeze(pns(Iak,:,:)),squeeze(mld));
        [pns(Iak,:,:)] = cut_off(squeeze(pns(Iak,:,:)),squeeze(pns(Iak,:,:)),squeeze(mld));
        
        display('optimizing density surface');
        tic
        [sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:),ithist] = optimize_surface_exact(s,ct,p,g,n2,sns(Iak,:,:),ctns(Iak,:,:),pns(Iak,:,:),e1t,e2t,lats,longs,settings);
        display(['optimizing density surface took ',int2str(toc),' seconds']);
        
        % save variables for postprocessing 
        save('ansu_hist.mat', 'ithist')
        
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
end
save (['omega_gk2.mat'])