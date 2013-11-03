clear all
close all
f1='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp229/data/k_zc.nc';
f2='/home/nfs/z3439823/eclipse/workspace/ansu/k_zc.nc';
dr1=ncread(f1,'k_zc');
dr2=ncread(f2,'k_zc');

max(abs(dr1(:)-dr2(:)))

diff=dr1'-dr2';
% h=imagesc(abs(diff))
% colorbar()
% set(gca,'YDir','normal')
% set(h,'alphadata',~isnan(diff))

f1='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp229/data/F.nc';
f2='/home/nfs/z3439823/eclipse/workspace/ansu/F.nc';
dr1=ncread(f1,'F');
dr2=ncread(f2,'F');

max(abs(dr1(:)-dr2(:)))

f1='/home/nfs/z3439823/eclipse/workspace/ansu/s.nc';
dr1=ncread(f1,'s');
f2='/home/nfs/z3439823/eclipse/workspace/ansu/s_.nc';
dr2=ncread(f2,'s_');

diff=dr1'-dr2';
h=imagesc(abs(diff))
colorbar()
set(gca,'YDir','normal')
set(h,'alphadata',~isnan(diff))
