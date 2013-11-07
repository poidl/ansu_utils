function n2ns=n2(s,ct,p,pns)

% calculate N^2 by taking S and Theta from 2 points upward and two points downward from the surface.
% take the mean pressure between these two points

[zi,yi,xi]=size(s);
s=s(1:zi,:);
ct=ct(1:zi,:);
p=p(1:zi,:);


[n2,pmid]=gsw_Nsquared(s,ct,p);

n2=p(1:end-1,:);

ii=bsxfun(@times,1:yi*xi,ones(zi-1,1));
pns_stacked=pns(ii); 

up=pmid<pns_stacked ;
kup=sum(up,1);
%notnan=~isnan(pmid(1,:,:));
kup(kup==0)=nan;
%kup(kup==0 & ~notnan)=nan;
%kup(kup==0 & notnan)=1; % for points between free surface and uppermost mid-pressure, set n2 equal to value at uppermost mid-pressure

kup=kup+(zi-1)*[0:xi*yi-1]; % 3d

bad=isnan(kup);
kup(bad)=2; % dummy, remove later

n2_1=n2(kup);
n2_2=n2(kup+1);
pmid_1=pmid(kup);
pmid_2=pmid(kup+1);

n2_1(bad)=nan; n2_2(bad)=nan;
pmid_1(bad)=nan; pmid_2(bad)=nan; 

dp=(pns(:)-pmid_1')./(pmid_2'-pmid_1');
n2ns=n2_1+(n2_2-n2_1).*dp';

n2ns=reshape(n2ns,[yi,xi]);

h=imagesc(n2ns)
set(gca,'YDir','normal')
set(h,'alphadata',~isnan(n2ns))
colorbar
title('N^2')
print('-dpdf','-r200',['figures/n.pdf'])

% edges=log(1:10:1800); 
% h=histc(n2ns(:),edges)
% bar(1:10:1800, h)
hist(n2ns(:),100)
title('N^2')
xlim(1e-3*[-1 1])
print('-dpdf','-r200',['figures/n2_hist.pdf'])
 
n2ns(n2ns<=0)=nan;
 
h=imagesc(sqrt(n2ns))
set(gca,'YDir','normal')
set(h,'alphadata',~isnan(n2ns))
colorbar
title('N')
print('-dpdf','-r200',['figures/n.pdf'])

end