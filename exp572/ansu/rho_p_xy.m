function [kx_ns,ky_ns]=rho_p_xy(s,ct,p,pns);


user_input;

[zi,yi,xi]=size(s);
s=s(1:zi,:);
ct=ct(1:zi,:);
p=p(1:zi,:);

ro=gsw_rho(s,ct,p);

rho_p=ro.*gsw_kappa(s,ct,p);

rho_p=reshape(rho_p,[zi,yi,xi]);
p=reshape(p,[zi,yi,xi]);

rho_px=circshift(rho_p, [0 0 -1])-rho_p;
px=0.5*(circshift(p, [0 0 -1])+p);
pnsx=0.5*(circshift(pns, [0 -1])+pns);
if ~zonally_periodic;
    rho_px(:,:,xi) = nan;
    px(:,:,xi) = nan;
    pnsx(:,xi) = nan;    
end
rho_py=circshift(rho_p, [0 -1 0])-rho_p;
py=0.5*(circshift(p, [0 -1 0])+p);
pnsy=0.5*(circshift(pns, [-1 0])+pns);
rho_py(:,yi,:) = nan;
py(:,yi,:) = nan;
pnsy(yi,:) = nan;

rho_px=rho_px(:,:);
rho_py=rho_py(:,:);
px=px(:,:);
py=py(:,:);

ii=bsxfun(@times,1:yi*xi,ones(zi,1));
pnsx_stacked=pnsx(ii); 
pnsy_stacked=pnsy(ii);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x
upx=px<pnsx_stacked;
kupx=sum(upx,1);
kupx(kupx==0)=nan;

% kup(kup==1)=nan; % recipe works only 2 points away from bottom and surface
% kup(kup==(zi-1))=nan;
kupx=kupx+zi*[0:xi*yi-1]; % 3d

bad=isnan(kupx);
kupx(bad)=2; % dummy, remove later

rho_px1=rho_px(kupx);
rho_px2=rho_px(kupx+1);
p1x=px(kupx);
p2x=px(kupx+1);

rho_px1(bad)=nan; rho_px2(bad)=nan;
p1x(bad)=nan; p2x(bad)=nan;

dpx=(pnsx(:)-p1x')./(p2x'-p1x');
kx_ns=rho_px1+(rho_px2-rho_px1).*dpx';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y
upy=py<pnsy_stacked;
kupy=sum(upy,1);
kupy(kupy==0)=nan;

% kup(kup==1)=nan; % recipe works only 2 points away from bottom and surface
% kup(kup==(zi-1))=nan;
kupy=kupy+zi*[0:xi*yi-1]; % 3d

bad=isnan(kupy);
kupy(bad)=2; % dummy, remove later

rho_py1=rho_py(kupy);
rho_py2=rho_py(kupy+1);
p1y=py(kupy);
p2y=py(kupy+1);

rho_py1(bad)=nan; rho_py2(bad)=nan;
p1y(bad)=nan; p2y(bad)=nan;

dpy=(pnsy(:)-p1y')./(p2y'-p1y');
ky_ns=rho_py1+(rho_py2-rho_py1).*dpy';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kx_ns=reshape(kx_ns,[yi,xi]);
ky_ns=reshape(ky_ns,[yi,xi]);

%  h=imagesc(sqrt(n2ns))
%  set(gca,'YDir','normal')
%  set(h,'alphadata',~isnan(n2ns))
%  colorbar
%  title('N')
%  print('-dpdf','-r200',['figures/n.pdf'])

end
