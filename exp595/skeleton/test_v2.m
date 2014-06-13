close all
clear all

load('data/bb_coords_surface.mat','ilat_','ilon_')
dlat=diff(ilat_);
dlat=dlat~=0;
dlat=[true,dlat];
dlon=diff(ilon_);
dlon=dlon~=0;
dlon=[true,dlon];
dl=dlat|dlon;
uilat_surf=ilat_(dl);
uilon_surf=ilon_(dl);
nbb_surf=sum(dl);

load('data/bb_coords_bottom.mat','ilat_','ilon_')
dlat=diff(ilat_);
dlat=dlat~=0;
dlat=[true,dlat];
dlon=diff(ilon_);
dlon=dlon~=0;
dlon=[true,dlon];
dl=dlat|dlon;
uilat_bot=ilat_(dl);
uilon_bot=ilon_(dl);
nbb_bot=sum(dl);

nk_surface=[];
nk_bottom=[];
ps=[];


for ii=1:nbb_bot
    s1=ncread(['data/nc/bottom/sns3d/sns3d',num2str(ii),'.nc'],'sns3d');
    ct1=ncread(['data/nc/bottom/ctns3d/ctns3d',num2str(ii),'.nc'],'ctns3d');
    p1=ncread(['data/nc/bottom/pns3d/pns3d',num2str(ii),'.nc'],'pns3d');
    s1=permute(s1,[3 2 1]);
    ct1=permute(ct1,[3 2 1]);
    p1=permute(p1,[3 2 1]); 

    if isempty(ps)
        ss=s1;
        cts=ct1;        
        ps=p1;
    else
        ss=cat(1,ss,s1);
        cts=cat(1,cts,ct1);
        ps=cat(1,ps,p1);
    end      
    
    nk_bottom=[nk_bottom,size(p1,1)];
end

for ii=nbb_bot+1:nbb_bot+nbb_surf
    
    s1=ncread(['data/nc/surface/sns3d/sns3d',num2str(ii),'.nc'],'sns3d');
    ct1=ncread(['data/nc/surface/ctns3d/ctns3d',num2str(ii),'.nc'],'ctns3d');
    p1=ncread(['data/nc/surface/pns3d/pns3d',num2str(ii),'.nc'],'pns3d');
    s1=permute(s1,[3 2 1]);
    ct1=permute(ct1,[3 2 1]);
    p1=permute(p1,[3 2 1]);  
    
    ss=cat(1,ss,s1);
    cts=cat(1,cts,ct1);
    ps=cat(1,ps,p1);  
    
    nk_surface=[nk_surface,size(p1,1)];
end

nk_surface=cumsum(nk_surface);
nk_bottom=cumsum(nk_bottom);

load('data/bb_coords_surface.mat')
ilon_s=ilon_;
ilat_s=ilat_;
load('data/bb_coords_bottom.mat')
ilon=[ilon_,ilon_s];
ilat=[ilat_,ilat_s];
load('data/p_bb_surface.mat')
p_bb_s=p_bb;
load('data/p_bb_bottom.mat')
pbb=[p_bb,p_bb_s];

uilon=[uilon_bot,uilon_surf];
uilat=[uilat_bot,uilat_surf];
nbb=length(uilon); % number of backbones

[ns,ny,nx]=size(ss);

ss_new=nan*ones(size(ss));
cts_new=nan*ones(size(ss));
ps_new=nan*ones(size(ss));
ilat_new=nan*ones(1,ns);
ilon_new=nan*ones(1,ns);

ilat_new_bb=nan*ones(1,ns);
ilon_new_bb=nan*ones(1,ns);

i1=1;
%keyboard
for ii=1:nbb
    ilo=uilon(ii);
    ila=uilat(ii);
    
    sdef=~isnan(ps(:,ila,ilo)); % true if surface intersects with backbone
    
    if any(sdef)
        ss_bb=ss(sdef,:,:);
        cts_bb=cts(sdef,:,:);
        ps_bb=ps(sdef,:,:);
        ilat_bb=ilat(sdef);
        ilon_bb=ilon(sdef);

        ps(sdef,:,:)=nan; % these surfaces are stored; delete

        [ilat_bb,ilon_bb,ss_bb,cts_bb,ps_bb]=sortit(ilat_bb,ilon_bb,ss_bb,cts_bb,ps_bb);

        i2=i1+sum(sdef)-1;
        %keyboard
        ss_new(i1:i2,:,:)=ss_bb;
        cts_new(i1:i2,:,:)=cts_bb;
        ps_new(i1:i2,:,:)=ps_bb;
        ilat_new(i1:i2)=ilat_bb;
        ilon_new(i1:i2)=ilon_bb;
        
        ilat_new_bb(i1:i2)=ila;
        ilon_new_bb(i1:i2)=ilo;

        i1=i2+1;
    end
end


dlat=diff(ilat_new_bb);
dlat=dlat~=0;
dlat=[true,dlat];
dlon=diff(ilon_new_bb);
dlon=dlon~=0;
dlon=[true,dlon];
dl=dlat|dlon;
uilat_surf=ilat_new_bb(dl);
uilon_surf=ilon_new_bb(dl);
nbb=sum(dl)

keyboard
sns3d=ss_new;
ctns3d=cts_new;
pns3d=ps_new;

% disp('final test: ')
% if checkit(pns3d)
%     error('not good')
% else
%     disp('success')
% end

save('data/bb_coords_final.mat','ilat','ilon')
save_netcdf(sns3d,'sns3d',['data/nc/sns3d_final.nc'])
save_netcdf(ctns3d,'ctns3d',['data/nc/ctns3d_final.nc'])
save_netcdf(pns3d,'pns3d',['data/nc/pns3d_final.nc'])

keyboard
