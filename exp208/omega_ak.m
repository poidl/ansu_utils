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

sa=s; clear s; % the _subs_ data saves sa in variable 's'

lats=lats(1,:,1); longs=squeeze(longs(1,1,:))';

ilon=find(longs>=225,1,'first');
ilat=find(lats>=-60,1,'first');

ctdash=zeros(size(ct));
ctdash(:,ilat,ilon)=linspace(1.0, 0.1, size(ct,1));
ct=ct+ctdash;

display('optimizing density surface');
tic

%dbstop in  optimize_surface_exact at 232 if ii==1220
%dbstop in optimize_surface_exact at 177


load('data/iteration_history_base.mat');
sns=squeeze(sns_hist(end,:,:));
ctns=squeeze(ctns_hist(end,:,:));
pns=squeeze(pns_hist(end,:,:));

[sns_i,ctns_i,pns_i ]= optimize_surface_exact(sa,ct,p,sns,ctns,pns);
display(['optimizing density surface took ',num2str(toc),' seconds']);

        
  
