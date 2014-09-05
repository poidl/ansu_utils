close all
clear all

sns3d=ncread('data/nc/sns3d_final.nc','sns3d');
ctns3d=ncread('data/nc/ctns3d_final.nc','ctns3d');
pns3d=ncread('data/nc/pns3d_final.nc','pns3d');

sns3d=permute(sns3d,[3 2 1]);
ctns3d=permute(ctns3d,[3 2 1]);
pns3d=permute(pns3d,[3 2 1]);
    
[ns,ny,nx]=size(pns3d);



cnt=0;
for kk=1:ns-1
    m  =~isnan(pns3d(kk,:,:)); % mask
    if sum(m(:))~=0;
        ex=kk;
        for ll=kk+1:ns
            m_b=~isnan(pns3d(ll,:,:)); % mask below
            if any(m(:) & m_b(:)); % common lateral grid points?
                m=m|m_b; % union
                ex=[ex,ll];
            end
        end

        cnt=cnt+1;

        ss=sns3d(ex,:,:);
        cts=ctns3d(ex,:,:);
        ps=pns3d(ex,:,:);

        pns3d(ex,:,:)=nan; % delete

        save_netcdf(ss,'ss',['data/nc/lateral/ss',num2str(cnt),'.nc'])
        save_netcdf(cts,'cts',['data/nc/lateral/cts',num2str(cnt),'.nc'])
        save_netcdf(ps,'ps',['data/nc/lateral/ps',num2str(cnt),'.nc'])
    end
end