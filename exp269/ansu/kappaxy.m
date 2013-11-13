function [kx_ns,ky_ns]=kappaxy(s,ct,p,pns);


user_input;

[zi,yi,xi]=size(s);
s=s(1:zi,:);
ct=ct(1:zi,:);
p=p(1:zi,:);


kappa=gsw_kappa(s,ct,p);

kappa=reshape(kappa,[zi,yi,xi]);
p=reshape(p,[zi,yi,xi]);

kappax=circshift(kappa, [0 0 -1])-kappa;
px=0.5*(circshift(p, [0 0 -1])+p);
pnsx=0.5*(circshift(pns, [0 -1])+pns);
if ~zonally_periodic;
    kappax(:,:,xi) = nan;
    px(:,:,xi) = nan;
    pnsx(:,xi) = nan;    
end
kappay=circshift(kappa, [0 -1 0])-kappa;
py=0.5*(circshift(p, [0 -1 0])+p);
pnsy=0.5*(circshift(pns, [-1 0])+pns);
kappay(:,yi,:) = nan;
py(:,yi,:) = nan;
pnsy(yi,:) = nan;

kappax=kappax(:,:);
kappay=kappay(:,:);
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

kappax1=kappax(kupx);
kappax2=kappax(kupx+1);
p1x=px(kupx);
p2x=px(kupx+1);

kappax1(bad)=nan; kappax2(bad)=nan;
p1x(bad)=nan; p2x(bad)=nan;

dpx=(pnsx(:)-p1x')./(p2x'-p1x');
kx_ns=kappax1+(kappax2-kappax1).*dpx';

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

kappay1=kappay(kupy);
kappay2=kappay(kupy+1);
p1y=py(kupy);
p2y=py(kupy+1);

kappay1(bad)=nan; kappay2(bad)=nan;
p1y(bad)=nan; p2y(bad)=nan;

dpy=(pnsy(:)-p1y')./(p2y'-p1y');
ky_ns=kappay1+(kappay2-kappay1).*dpy';

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
