% subsample and convert sp to sa
clear all

% Load selected variables.
vars = {'s','ct','p','lats','longs'};
load('/home/z3439823/mymatlab/omega/data_paul/gk_ak_gamma.mat', vars{:})

sk=1;
xskip=sk;
yskip=sk;
zskip=sk;

% subsample
s=s(1:zskip:end,1:yskip:end,1:xskip:end);
ct=ct(1:zskip:end,1:yskip:end,1:xskip:end);
p=p(1:zskip:end,1:yskip:end,1:xskip:end);
lat=lats(1:yskip:end);
lon=longs(1:xskip:end);

ss=size(s);
lon=repmat(permute(lon,[3,2,1]),[ss(1),ss(2),1]);
lat=repmat(permute(lat,[3,1,2]),[ss(1),1,ss(3)]);
    
% convert sp to sa
s=gsw_SA_from_SP(s,p,lon,lat);

lats=lat;
longs=lon;

save('data/input_data.mat',vars{:})

