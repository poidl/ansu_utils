clear all;
close all;

runs=344:442;
%runs=344:346;

addpath(genpath('../exp344/ansu')) % for find_regions

fname=['../exp344/data/stationindex.mat'];
load(fname);
disp(istation)

fname=['../exp344/data/iteration_history.mat'];
load(fname,'pns_hist','sns_hist','ctns_hist');

nit=size(pns_hist,1);
[yi,xi]=size(squeeze(pns_hist(1,:,:)));

sns_3d=nan*ones(length(runs),yi,xi);
ctns_3d=nan*ones(length(runs),yi,xi);
pns_3d=nan*ones(length(runs),yi,xi);

p_istation=nan*ones(1,length(runs));
ct_istation=nan*ones(1,length(runs));
s_istation=nan*ones(1,length(runs));

if 1
    cnt=0;
    for ii=runs
        if mod(ii,5)==0
            disp(['reading exp',num2str(ii)])
        end
        cnt=cnt+1;
        fname=['../exp',num2str(ii),'/data/iteration_history.mat'];
        load(fname,'pns_hist');
        load(fname,'ctns_hist');
        load(fname,'sns_hist');
        pns=squeeze(pns_hist(nit,:,:));
        sns=squeeze(sns_hist(nit,:,:));
        ctns=squeeze(ctns_hist(nit,:,:));
        setnan=true(size(sns));
        regions=find_regions(sns);
        found=false;
        for iregion=1:length(regions)
            region=regions{iregion};
            if ismember(istation,region)
                setnan(region)=false;
                pns(setnan)=nan;
                sns(setnan)=nan;
                ctns(setnan)=nan;
                pns_3d(cnt,:,:)=pns;
                sns_3d(cnt,:,:)=sns;
                ctns_3d(cnt,:,:)=ctns;
                p_istation(cnt)=pns(istation);
                ct_istation(cnt)=ctns(istation);
                s_istation(cnt)=sns(istation);
                found=true;
            end
        end
        if found==false
            error('hoitaus')
        end
    end
    save skeleton.mat p_istation ct_istation s_istation pns_3d sns_3d ctns_3d
else
    load skeleton.mat
    for ii=1:nit-1
        both=~isnan(surfvar(nit,:)) & ~isnan(surfvar(nit,:));
        if any( surfvar(nit,both)>surfvar(nit+1,both))
            disp(['error in exp',num2str(runs(ii))])
        end
    end
end





% 
% vp=squeeze(surfvar(1,:,:));
% h=imagesc( vp )
% set(h,'alphadata',~isnan(vp)) % white nans
% set(gca,'YDir','normal')
% colorbar()

