clear all
close all

load('data/bb_coords_final.mat')

%ncread('data/nc/sns3d_final.nc','sns3d')
%ncread('data/nc/ctns3d_final.nc','ctns3d')
pns3d=ncread('data/nc/pns3d_final.nc','pns3d');

[ns,ny,nx]=size(pns3d);

for kk=1:ns
    cnt=0;
    for ii=1:length(ilon)
        if ~isnan(pns3d(kk,ilat(ii),ilon(ii)))
            cnt=cnt+1;
        end
    end
    disp(['level ',num2str(kk),': ',num2str(cnt)])
end

