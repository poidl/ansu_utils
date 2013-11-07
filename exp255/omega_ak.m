clear all
close all

params; % load some parameters

% load('data/input_data.mat')
% for ii=1:size(s,3)
%     for jj=1:size(s,2)
%         kk=find(isnan(s(:,jj,ii)),1,'first');
%         s(kk:end,jj,ii)=nan;
%         ct(kk:end,jj,ii)=nan;
% %        p(kk:end,jj,ii)=nan;
%     end
% end
load('data/idealized_01.mat')

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
    
    disp('not relevant here')
    
else
    % for data set idealized_01.mat (constant salinity), isothermals are
    % pot. dens. surfaces and neutral surfaces
    [zi,yi,xi]=size(s);
%     isouth=find(ct(:,1,1)==1,1,'first');
%     inorth=find(ct(:,2,1)==1,1,'first');
    isouth=41;
    inorth=61;
    
    ctns=nan*ones([yi,xi]);
    pns=nan*ones([yi,xi]);
    
    ctns(1,:)=ct(isouth,1,1);
    ctns(2,:)=ct(inorth+1,2,1);
    pns(1,:)=p(isouth,1,1);
    pns(2,:)=p(inorth+1,2,1);
    
    sns=35+0*ctns;
    
    ctns(:,1)
    pns(:,1)
%     ptarget=p(ip,1,1)
%     sns=squeeze( s(ip,:,:));
%     ctns=squeeze( ct(ip,:,:));
%     pns=squeeze( p(ip,:,:));
%     pns(isnan(sns))=nan;
%     if 0;
%         ipert=-20;
%         sns(2,:)=squeeze( s(ip+ipert,2,:));
%         ctns(2,:)=squeeze( ct(ip+ipert,2,:));
%         pns(2,:)=squeeze( p(ip+ipert,1,1));          
%         pns(isnan(sns))=nan;
%     end
    
end



display('optimizing density surface');
tic

%dbstop in  optimize_surface_exact at 232 if ii==1220
%dbstop in optimize_surface_exact at 333

%save_netcdf(sa,ct,p,sns,ctns,pns);
[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:)] = optimize_surface_exact(s,ct,p,sns,ctns,pns);
display(['optimizing density surface took ',num2str(toc),' seconds']);

        
  
