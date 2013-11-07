function n2ns=n2(s,ct,p,pns)

% calculate N^2 by taking S and Theta from 2 points upward and two points downward from the surface.
% take the mean pressure between these two points

[zi,yi,xi]=size(s);
s=s(1:zi,:);
ct=ct(1:zi,:);
p=p(1:zi,:);
ii=bsxfun(@times,1:yi*xi,ones(zi,1));
pns_stacked=pns(ii); 

up=p<pns_stacked ;
kup=sum(up,1);
kup(kup==0)=nan;

kup(kup==1)=nan; % recipe works only 2 points away from bottom and surface
kup(kup==(zi-1))=nan;
kup=kup+zi*[0:xi*yi-1]; % 3d

bad=isnan(kup);
kup(bad)=2; % dummy, remove later

s1=s(kup-1);
s2=s(kup+2);
ct1=ct(kup-1);
ct2=ct(kup+2);
p1=p(kup-1);
p2=p(kup+2);

s1(bad)=nan; s2(bad)=nan;
ct1(bad)=nan; ct2(bad)=nan;
p1(bad)=nan; p2(bad)=nan;

smid=0.5*(s1+s2);
ctmid=0.5*(ct1+ct2);
pmid=0.5*(p1+p2);

dp=(p2-p1);
rz=(1./dp).*(gsw_rho(s1,ct1,pmid)-gsw_rho(s2,ct2,pmid));

r=gsw_rho(smid,ctmid,pmid);

rz=reshape(rz,[yi,xi]);
r=reshape(r,[yi,xi]);
n2ns=-9.81*(rz./r);

%  h=imagesc(sqrt(n2ns))
%  set(gca,'YDir','normal')
%  set(h,'alphadata',~isnan(n2ns))
%  colorbar
%  title('N')
%  print('-dpdf','-r200',['figures/n.pdf'])

end