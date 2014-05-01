clear all
close all

params; % load some parameters

load('data/input_data.mat')
% for ii=1:size(s,3)
%     for jj=1:size(s,2)
%         kk=find(isnan(s(:,jj,ii)),1,'first');
%         if ~isempty(kk)
%             if kk(1)==1 & 
%                 disp('ice?')
%             end
%         end
%         s(kk:end,jj,ii)=nan;
%         ct(kk:end,jj,ii)=nan;
% %        p(kk:end,jj,ii)=nan;
%     end
% end

% load('data/b.mat')
% b=b(1:4:end,:,:);
% for ii=1:size(s,3)
%     for jj=1:size(s,2)
%         if all(isnan(b(:,jj,ii)))
%             s(:,jj,ii)=nan;
%             ct(:,jj,ii)=nan;
%         end
%     end
% end


% the lats/longs are only used to calculate epsilon in diagnose_and_write()
% they are not necessary to calculate omega surfaces
lats=squeeze(lats(1,:,:)); longs=squeeze(longs(1,:,:));
[dy,dx]=scale_fac(lats,longs);
save dxdy.mat dx dy

vsurf=var_on_surf_stef(va,p,surf)

lat=lats(:,1);
lon=longs(1,:);
[mini,ilat]=min(abs(lat+16)); % Jackett & McDougall 97: 16 South 188 North
[mini,ilon]=min(abs(lon-188));
[zi,yi,xi]=size(s);
istation=ilat+yi*(ilon-1);
save('data/stationindex.mat','istation')

sa=s; clear s; % the _subs_ data saves sa in variable 's'

% load('../exp266/data/iteration_history.mat','sns_hist');
% regions=find_regions(squeeze(sns_hist(end,:,:)));
% 
% [zi,yi,xi]=size(sa);
% imaxregion=0; % index of largest region
% npmaxreg=0; % number of points in largest region
% for i=1:length(regions)
%     if length(regions{i})>npmaxreg;
%         npmaxreg=length(regions{i});
%         imaxreg=i;
%     end
% end
% setnan=true(1,xi*yi);
% setnan(regions{imaxreg})=false;
% sa(:,setnan)=nan;
% ct(:,setnan)=nan;
% %p(:,setnan)=nan;


[zi,yi,xi] = size(sa);

sns_i = nan(length(nlevels),yi,xi);
ctns_i = nan(length(nlevels),yi,xi);
pns_i = nan(length(nlevels),yi,xi);

Iak = 1;

if initial_surface_at_constant_pressure;
    ip=find(p(:,1,1)>=initial_pressure,1,'first');
    ptarget=p(ip,1,1)
    sns=squeeze(sa(ip,:,:));
    ctns=squeeze(ct(ip,:,:));
    pns=squeeze(p(ip,:,:));
    pns(isnan(sns))=nan;
    
else
    glevels_Iak=1e3+glevels(Iak);
    [zi,yi,xi]=size(sa);
    
    s_2=sa(:,:);
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
        %keyboard
        F=gsw_rho(s_2,ct_2,p_r*ones(size(s_2)))-glevels_Iak;
        %dbstop in root_core at 29
        [s_2,ct_2,p_2,sns_out,ctns_out,pns_out, inds, fr]=root_core(F,inds,refine_ints,s_2,ct_2,p_2,sns_out,ctns_out,pns_out);
        
        if all(~fr) % break out of loop if all roots have been found
            break
        end
        
    end
    sns=sns_out;
    ctns=ctns_out;
    pns=pns_out;
end


display('optimizing density surface');
tic
%dbstop in optimize_surface_exact at 624
%dbstop in optimize_surface_exact at 225
%dbstop in  optimize_surface_exact at 232 if ii==1220
%dbstop in optimize_surface_exact at 177
%dbstop in optimize_surface_exact at 104
%dbstop in optimize_surface_exact at 328

load('data/mask.mat')
sns(~mask)=nan;
ctns(~mask)=nan;
pns(~mask)=nan;

load('data/stationindex.mat');
setnan=true(size(sns));
regions=find_regions(sns);
for iregion=1:length(regions)
    region=regions{iregion};
    if ismember(istation,region)
        setnan(region)=false;        
        pns(setnan)=nan;
        sns(setnan)=nan;
        ctns(setnan)=nan;      
    end
end

[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:)] = optimize_surface_exact(sa,ct,p,sns,ctns,pns);
display(['optimizing density surface took ',num2str(toc),' seconds']);

        
        
  
