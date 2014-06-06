lat=lat(:,1);
lon=lon(1,:);
[mini,ilat]=min(abs(lat+16)); % Jackett & McDougall 97: 16 South 188 East
[mini,ilon]=min(abs(lon-188));
[zi,yi,xi]=size(sa);
istation=ilat+yi*(ilon-1);
save('data/stationindex.mat','istation')

% artificial bottom for testing
sa(end,:,:)=nan;
ct(end,:,:)=nan;
bot=squeeze(sum(~isnan(sa),1)+1);
bot3d=bot(:)'+zi*(0:yi*xi-1);
bottom=p(bot3d);
bottom=reshape(bottom,[yi xi]);
land=isnan(sa(1,:));
bottom(land)=nan;

save_netcdf(bottom,'bottom','data/bottom.nc')

depth_bb=p(bot3d(ilat+yi*(ilon-1))-1); % bottom depth at backbone
ns=3;
nstmp=ns+2;
pbb=linspace(0,depth_bb,nstmp); 
pbb=pbb(2:end-1);


sns3d=nan*ones(ns,yi,xi);
ctns3d=sns3d;
pns3d=sns3d;

if 0
    for kk=1:ns

        point=[pbb(kk) ilat ilon];
        display(['optimizing pressure ',num2str(pbb(kk)),' on backbone (depth: ',num2str(depth_bb),')']);
        [sns3d(kk,:,:),ctns3d(kk,:,:),pns3d(kk,:,:)] =optimize_surface_at_point(sa,ct,p,point);

    end
    save('data/omega_3d.mat','sns3d','ctns3d','pns3d')
else
    load('data/omega_3d.mat')
end


dp=diff(pbb); dp=dp(1);

p_1=squeeze(pns3d(1,:,:));
p_2=squeeze(pns3d(2,:,:));

p_1(~isnan(p_2))=nan;
dist=bottom-p_1;
todo=dist>dp;
save_netcdf(double(todo),'todo','data/todo.nc')
keyboard


p_lo=squeeze(pns3d(end,:,:));
bottom-p_lo;



