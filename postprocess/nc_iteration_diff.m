clear all
close all
f1='/home/nfs/z3439823/eclipse/workspace/ansu/pns2.nc';
f2='/home/nfs/z3439823/eclipse/workspace/ansu/pns3.nc';
dr1=ncread(f1,'pns');
dr2=ncread(f2,'pns');

max(abs(dr1(:)-dr2(:)))

diff=dr1'-dr2';
h=imagesc(diff)
colorbar()
set(gca,'YDir','normal')

set(h,'alphadata',~isnan(diff))
cmp=colormap(hot(128));
cmp2=cmp(find(cmp(:,1)==1,1,'first')-15:end,:);
colormap([fliplr(cmp2);flipud(cmp2)]) ;
maxi=max(abs(diff(:)));
caxis([-maxi maxi])

