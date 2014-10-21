function [n2ns,n2nsx,n2nsy]=get_n2ns(pns,s,ct,p)


[ny,nx]=size(pns);
[n2,~]=n2_smooth(s,ct,p);
n2(n2(:)<=1e-6)=1e-6;

pns_zonal=nanmean(pns,2);
pns_zonal=repmat(pns_zonal,[1 nx]);
n2ns=var_on_surf_stef(n2,p,pns_zonal);

n2nsx=0.5*(n2ns+circshift(n2ns,[0 -1]));
n2nsx(:,end)=n2nsx(:,end-1); % sloppy 

n2nsy=0.5*(n2ns+circshift(n2ns,[-1 0]));
n2nsy(end,:)=n2nsy(end-1,:); % sloppy   

