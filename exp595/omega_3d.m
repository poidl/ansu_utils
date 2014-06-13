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
p_max=zeros(yi,xi);
p_min=last;
p_max(land)=nan;
dist=last-p_max;

[md,imax]=max(dist(:));

p1=0;
p2=p(last3d(ilat+yi*(ilon-1))-1); % depth of deepest data point at backbone
dp=500;
pbb=(dp:dp:p2);

ilat_=[];
ilon_=[];
p_bb=[];

cnt=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill out bottom holes
while 1;
    disp(['backbone Nr. ',num2str(cnt)])
    [sns3d,ctns3d,pns3d] = ribs(ilat,ilon,p1,p2,pbb,sa,ct,p);
    
    ns=size(sns3d,1);
    for kk=1:ns
        ps=squeeze(pns3d(kk,:,:));
        deeper=(ps>p_max);
        shallower=(ps<p_min);
        p_max(deeper)=ps(deeper);
        p_min(shallower)=ps(shallower);
    end

    dist=last-p_max;
    
    ilat_=[ilat_,ilat*ones(1,ns)];
    ilon_=[ilon_,ilon*ones(1,ns)];
    p_bb=[p_bb,repmat([p1;p2],[1 ns])];
    save_netcdf(p_max,'p_max',['data/nc/bottom/p_max/p_max',num2str(cnt),'.nc'])
    save_netcdf(dist,'dist',['data/nc/bottom/dist/dist',num2str(cnt),'.nc'])
    
    save_netcdf(sns3d,'sns3d',['data/nc/bottom/sns3d/sns3d',num2str(cnt),'.nc'])
    save_netcdf(ctns3d,'ctns3d',['data/nc/bottom/ctns3d/ctns3d',num2str(cnt),'.nc'])
    save_netcdf(pns3d,'pns3d',['data/nc/bottom/pns3d/pns3d',num2str(cnt),'.nc'])
    
    [md,imax]=max(dist(:));
    if md>dp
        [ilat,ilon]=ind2sub(size(dist),imax);
        istation=ilat+yi*(ilon-1);
        save('data/stationindex.mat','istation')

        p1=p_max(ilat,ilon);
        p2=last(ilat,ilon); % depth of deepest data point at backbone
        pbb=(p1+dp:dp:p2);
        cnt=cnt+1;
    else
        break
    end
end
save('data/bb_coords_bottom.mat','ilat_','ilon_')
save('data/p_bb_bottom.mat','p_bb')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill out surface holes

dist=p_min;
[md,imax]=max(dist(:));

ilat_=[];
ilon_=[];
p_bb=[];
while md>dp;
    cnt=cnt+1;
    disp(['backbone Nr. ',num2str(cnt)])    
    
    [ilat,ilon]=ind2sub(size(dist),imax);
    istation=ilat+yi*(ilon-1);
    save('data/stationindex.mat','istation')

    p1=0;
    p2=p_min(ilat,ilon);
    pbb=(p1+dp:dp:p2);
    
    [sns3d,ctns3d,pns3d] = ribs(ilat,ilon,p1,p2,pbb,sa,ct,p);
    
    ns=size(sns3d,1);
    for kk=1:ns
        ps=squeeze(pns3d(kk,:,:));
        shallower=(ps<p_min);
        p_min(shallower)=ps(shallower);
    end
    
    ilat_=[ilat_,ilat*ones(1,ns)];
    ilon_=[ilon_,ilon*ones(1,ns)];
    p_bb=[p_bb,repmat([p1;p2],[1 ns])];
    save_netcdf(p_min,'p_min',['data/nc/surface/p_min/p_min',num2str(cnt),'.nc'])
    save_netcdf(dist,'dist',['data/nc/surface/dist/dist',num2str(cnt),'.nc'])
    
    save_netcdf(sns3d,'sns3d',['data/nc/surface/sns3d/sns3d',num2str(cnt),'.nc'])
    save_netcdf(ctns3d,'ctns3d',['data/nc/surface/ctns3d/ctns3d',num2str(cnt),'.nc'])
    save_netcdf(pns3d,'pns3d',['data/nc/surface/pns3d/pns3d',num2str(cnt),'.nc'])
    
    dist=p_min;
    [md,imax]=max(dist(:));
end
save('data/bb_coords_surface.mat','ilat_','ilon_')
save('data/p_bb_surface.mat','p_bb')

keyboard


