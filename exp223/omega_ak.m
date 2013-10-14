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

[longs,lats] = meshgrid(longs,lats);

[zi,yi,xi] = size(s);


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
    sns_Iak=s(ip,:,:);
    ctns_Iak=ct(ip,:,:);
    pns_Iak=p(ip,:,:);
    pns_Iak(isnan(sns_Iak))=nan;
else
glevels_Iak=1e3+glevels(Iak);
delta=1e-11;
[zi,yi,xi]=size(s);

s_2=s(:,:);
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

        [final,fr,k_zc]=root_core(F,delta,stack);

        k_zc_3d=k_zc+stack*[0:size(F,2)-1]; % indices of flattened 3d-array where root has been found   

        sns_out(inds(final))=s_2(k_zc_3d(final)); % adjust surface where root has already been found
        ctns_out(inds(final)) =ct_2(k_zc_3d(final));
        pns_out(inds(final)) =p_2(k_zc_3d(final));
        inds=inds(fr); % points where surface has not been corrected

        if all(~fr) % break out of loop if all roots have been found
            break
        end

        k=k_zc_3d(fr);  % indices of flattened 3d-array where vertical resolution must be increased

        %keyboard
        ds_ =  ( s_2(k+1) - s_2(k))/refine_ints; % increase resolution in the vertical
        dt_ = (ct_2(k+1) - ct_2(k))/refine_ints;
        dp_ =  (p_2(k+1) - p_2(k))/refine_ints;

        ds_ =bsxfun(@times, ds_, [0:refine_ints]');
        dt_ = bsxfun(@times, dt_, [0:refine_ints]');
        dp_ = bsxfun(@times, dp_, [0:refine_ints]');

        s_2 =  bsxfun(@plus,s_2(k),ds_);
        ct_2 =  bsxfun(@plus,ct_2(k),dt_);
        p_2 =  bsxfun(@plus,p_2(k),dp_);

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
%dbstop in optimize_surface_exact at 229

load('../exp203/data/iteration_history.mat');
sns2=squeeze(sns_hist(2,:,:));
ctns2=squeeze(ctns_hist(2,:,:));
pns2=squeeze(pns_hist(2,:,:));

sns(isnan(sns2))=nan;
ctns(isnan(ctns2))=nan;
pns(isnan(pns2))=nan;

[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:)] = optimize_surface_exact(s,ct,p,sns,ctns,pns);
display(['optimizing density surface took ',num2str(toc),' seconds']);

        
  
