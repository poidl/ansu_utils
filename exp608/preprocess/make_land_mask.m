close all

[nz,ny,nx]=size(sa);
[ocean, n] = gamma_ocean_and_n(sa,ct,p,lon,lat);

%north-south
test1=ocean.*circshift(ocean,[-1 0]);
test2=ocean.*circshift(ocean,[0 -1]);

n= (test1==5 | test2==5);
setnan=n(:);

for kk=1:nz
    sa(kk,setnan)=nan;
    ct(kk,setnan)=nan;
end


% [ocean2, n2] = gamma_ocean_and_n(sa,ct,p,lon,lat);
% h=imagesc(ocean)
% set(h,'alphadata',~isnan(ocean))
% figure()
% h=imagesc(ocean2)
% set(h,'alphadata',~isnan(ocean2))
