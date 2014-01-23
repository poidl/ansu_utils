close all
clear all

load b_original.mat
lat1=lat;
lon1=lon; 
p1=p;
b1=b;
clear lat lon p b

% lat1=permute(lat1,[3 2 1]);
% lon1=permute(lon1,[3 2 1]);
% p1=permute(p1,[3 2 1]);
% b1=permute(b1,[3 2 1]);

load('/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp324/data/input_data.mat')

%b=interp3(p1,lat1,lon1,b1,p,lats,longs);
b=interp3(lat1,p1,lon1,b1,lats,p,longs);

save b.mat b
