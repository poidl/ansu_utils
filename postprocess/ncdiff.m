clear all
close all
f1='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp229/data/drho.nc';
f2='/home/nfs/z3439823/eclipse/workspace/ansu/drho.nc';
dr1=ncread(f1,'drho');
dr2=ncread(f2,'drho');

max(abs(dr1(:)-dr2(:)))

h=imagesc(abs(dr1'-dr2'))
colorbar()
set(gca,'YDir','normal')
set(h,'alphadata',~isnan(diff))