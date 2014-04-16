clear all
close all
% f1='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp229/data/sns(neighbour).nc';
% f2='/home/nfs/z3439823/eclipse/workspace/ansu/sns_(inbr).nc';
% dr1=ncread(f1,'matl');
% dr2=ncread(f2,'fort');
% 
% [dr1,i1]=sort(dr1);
% [dr2,i2]=sort(dr2);
% 
% plot(dr2,'*r')
% hold on
% plot(dr1,'*')

% f1='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp229/data/s(:,en).nc';
% f2='/home/nfs/z3439823/eclipse/workspace/ansu/s_(iwo,:).nc';%
%f1='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp229/data/bottle.nc';
%f2='/home/nfs/z3439823/eclipse/workspace/ansu/bottle.nc';
%f2='/home/nfs/z3439823/fortran/debug/bottle.nc';
 f1='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp229/data/cast.nc';
% f2='/home/nfs/z3439823/eclipse/workspace/ansu/cast.nc';
f2='/home/nfs/z3439823/fortran/debug/cast.nc';
 f1='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp229/data/Fneighbour.nc';
% f2='/home/nfs/z3439823/eclipse/workspace/ansu/Fneighbour.nc';
f1='/home/nfs/z3439823/eclipse/workspace/ansu/pns_dumb.nc';
f2='/home/nfs/z3439823/eclipse/workspace/ansu/pns.nc';


dr1=ncread(f1,'pns');
dr2=ncread(f2,'pns');

v1=dr1(:,1);
v2=dr2(:,1);
[trash,i1]=sort(v1);
[trash,i2]=sort(v2);


dr1=dr1(i1,:);
dr2=dr2(i2,:);

dr1(dr1<-0.1)=-0.1;
dr1(dr1>0.1)=0.1;
close all
h=imagesc(dr1)
set(h,'alphadata',~isnan(dr1)) % white nans
set(gca,'YDir','normal')
title('matl')
colorbar
caxis([-0.1 0.1])

dr2(dr2<-0.1)=-0.1;
dr2(dr2>0.1)=0.1;
figure()
h=imagesc(dr2)
set(h,'alphadata',~isnan(dr2)) % white nans
set(gca,'YDir','normal')
title('fort')
colorbar
caxis([-0.1 0.1])


