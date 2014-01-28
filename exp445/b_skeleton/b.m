clear all;
close all;

fname=['../exp344/data/input_data.mat'];
load(fname);

p=p(:,:);
s=s(:,:);
ct=ct(:,:);

load omega.mat 

omega=omega(:,:);

pmid=0.5*(p(1:end-1,:)+p(2:end,:));
r1=gsw_rho(s(1:end-1,:),ct(1:end-1,:),pmid);
r2=gsw_rho(s(2:end,:),ct(2:end,:),pmid);

dp=p(1:end-1,:)-p(2:end,:);

mid_rgN2=(r1-r2)./dp;

om_z=(omega(1:end-1,:)-omega(2:end,:))./dp;

b=om_z./mid_rgN2;

b=reshape(b,[zi,yi,xi]);

save b.mat b

delete b.nc        
nccreate('b.nc','b',...
         'Dimensions', {'k', zi, 'j', yi, 'i', xi})
ncwrite('b.nc','b', b);
