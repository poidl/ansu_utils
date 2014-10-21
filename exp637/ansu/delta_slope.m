function [sx,sy]=delta_slope(sns,ctns,pns,s,ct,p)

% north 
sn=circshift(s,[0 -1 0]);
ctn=circshift(ct,[0 -1 0]);
pn=circshift(p,[0 -1 0]);    
% south 
ss=circshift(s,[0 1 0]);
cts=circshift(ct,[0 1 0]);
ps=circshift(p,[0 1 0]); 

% east
se=circshift(s,[0 0 -1]);
cte=circshift(ct,[0 0 -1]);
pe=circshift(p,[0 0 -1]);    
% west
sw=circshift(s,[0 0  1]);
ctw=circshift(ct,[0 0 1]);
pw=circshift(p,[0 0 1]);

% height of neutral tangent at eastern (...) neighbour, gridded onto point
% from which tangent is emanating
[tr,tr,pe]=depth_ntp_simple(sns(:)',ctns(:)',pns(:)',se(:,:),cte(:,:),pe(:,:)); 
[tr,tr,pw]=depth_ntp_simple(sns(:)',ctns(:)',pns(:)',sw(:,:),ctw(:,:),pw(:,:));
[tr,tr,pn]=depth_ntp_simple(sns(:)',ctns(:)',pns(:)',sn(:,:),ctn(:,:),pn(:,:));
[tr,tr,ps]=depth_ntp_simple(sns(:)',ctns(:)',pns(:)',ss(:,:),cts(:,:),ps(:,:));

[yn,xn]=size(sns);
pe=reshape(pe,[yn,xn]);
pw=reshape(pw,[yn,xn]);
pn=reshape(pn,[yn,xn]);
ps=reshape(ps,[yn,xn]);
%keyboard
% dz east of grid point (gridded onto original point)
pnse=circshift(pns,[0,-1]);
dz_east=(pe-pnse);
% dz at grid point
pw=circshift(pw,[0 -1]); % grid onto point to which tangent is laid
dz_orig=(pw-pns); 
sx=0.5*(dz_east-dz_orig);

% dz north of grid point (gridded onto original point)
pnsn=circshift(pns,[-1,0]);
dz_north=(pn-pnsn);
% dz at grid point
ps=circshift(ps,[-1 0]); % grid onto point to which tangent is laid
dz_orig=(ps-pns); 
sy=0.5*(dz_north-dz_orig);



