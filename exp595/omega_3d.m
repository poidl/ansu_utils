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

last3d=bot3d-1; % deepest data point
last3d(last3d==0)=1; % dummy
last=p(last3d);
last=reshape(last,[yi xi]);
last(land)=nan;

save_netcdf(bottom,'bottom','data/nc/bottom.nc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_e=zeros(yi,xi);
p_e(land)=nan;
dist=last-p_e;

[md,imax]=max(dist(:));

p1=0;
p2=p(last3d(ilat+yi*(ilon-1))-1); % depth of deepest data point at backbone
dp=500;
pbb=(dp:dp:p2);
md=p2;

ilat_=[];
ilon_=[];

cnt=1;
while md>dp;
    disp(['backbone Nr. ',num2str(cnt)])
    [sns3d,ctns3d,pns3d] = ribs(ilat,ilon,p1,p2,pbb,sa,ct,p);
    
    ns=size(sns3d,1);
    for kk=1:ns
        ps=squeeze(pns3d(kk,:,:));
        deeper=(ps>p_e);
        p_e(deeper)=ps(deeper);
    end

    if cnt==80
        keyboard
    end
    dist=last-p_e;
    
    ilat_=[ilat_,ilat];
    ilon_=[ilon_,ilon];
    save_netcdf(p_e,'p_e',['data/nc/p_e/p_e',num2str(cnt),'.nc'])
    save_netcdf(dist,'dist',['data/nc/dist/dist',num2str(cnt),'.nc'])
    
    save_netcdf(sns3d,'sns3d',['data/nc/sns3d/sns3d',num2str(cnt),'.nc'])
    save_netcdf(ctns3d,'ctns3d',['data/nc/ctns3d/ctns3d',num2str(cnt),'.nc'])
    save_netcdf(pns3d,'pns3d',['data/nc/pns3d/pns3d',num2str(cnt),'.nc'])
    
    [md,imax]=max(dist(:));
    [ilat,ilon]=ind2sub(size(dist),imax);
    istation=ilat+yi*(ilon-1);
    save('data/stationindex.mat','istation')
    
    p1=p_e(ilat,ilon);
    p2=last(ilat,ilon); % depth of deepest data point at backbone
    pbb=(p1+dp:dp:p2);
    cnt=cnt+1;
end
save('data/bb_coords.mat','ilat_','ilon_')
keyboard


