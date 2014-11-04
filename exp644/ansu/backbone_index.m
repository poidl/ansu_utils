function [ibb_horz,i,j]=backbone_index(lon2d,lat2d)


lvec=[lon2d(:)-188,lat2d(:)+16]; 
[~,ibb_horz]=min(lvec(:,1).^2+lvec(:,2).^2);
[i,j]=ind2sub(size(lat2d),ibb_horz);

