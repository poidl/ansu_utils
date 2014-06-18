close all
clear all

load('data/bb_coords.mat','ilat_','ilon_')
dlat=diff(ilat_);
dlat=dlat~=0;
dlat=[true,dlat];
dlon=diff(ilon_);
dlon=dlon~=0;
dlon=[true,dlon];
dl=dlat|dlon;
nbb=sum(dl); % number of backbones
lbb=cumsum(dl); % backbone label
uilat=ilat_(dl); % lats/lons of unique backbones
uilon=ilon_(dl);

ilat=ilat_;
ilon=ilon_;


ps=[];

for ii=1:nbb
    s1=ncread(['data/nc/stage1/sns3d/sns3d',num2str(ii),'.nc'],'sns3d');
    ct1=ncread(['data/nc/stage1/ctns3d/ctns3d',num2str(ii),'.nc'],'ctns3d');
    p1=ncread(['data/nc/stage1/pns3d/pns3d',num2str(ii),'.nc'],'pns3d');
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
    
end

load('data/p_bb.mat')
up_bb=p_bb(:,dl);

[ns,ny,nx]=size(ss);

ss_new=nan*ones(size(ss));
cts_new=nan*ones(size(ss));
ps_new=nan*ones(size(ss));
ilat_new=nan*ones(1,ns);
ilon_new=nan*ones(1,ns);

ilat_new_bb=nan*ones(1,ns);
ilon_new_bb=nan*ones(1,ns);

dp=500;

i1=1;
%keyboard
for ii=1:nbb
    ilo=uilon(ii);
    ila=uilat(ii);
    p1=up_bb(1,ii); 
    p2=up_bb(2,ii);
    
    sdef=~isnan(ps(:,ila,ilo)); % true if surface intersects with backbone
    
%    if any(sdef)
    ss_bb=ss(sdef,:,:);
    cts_bb=cts(sdef,:,:);
    ps_bb=ps(sdef,:,:);

    is=find(sdef);      

    for jj=1:sum(sdef)
        js=is(jj);
        if lbb(js)~=ii

            jsa=is(jj-1);
            ilo_=ilon(js);
            ila_=ilat(js);
            keyboard
            dp_rms=sqrt(nanmean((ps(js,:)-ps(jsa,:)).^2));
        end
    end
%             ilo_=ilon(ii);
%             ila_=ilat(ii);
%             p1_=p_bb(1,ii);
%             p2_=p_bb(2,ii);
%             
%             ibb= lbb==lbb(ii);
%             pbb=ps(ibb,ilat_,ilo_);
%             
%             d1=abs(ps(ii,ila_,ilo_)-pbb(2));
%             d2=abs(ps(ii,ila_,ilo_)-p2_);
%             ep=1e-5;
%             
%             ibb= lbb==lbb(ii);
%             
%             pbb=ps(ibb,ilat_,ilo_);
%             if d1<ep
%                 keyboard
%                 p_bb(1,ii)=pbb(2,ii);
%             elseif d2<ep
%                 keyboard
%                 p_bb(2,ii)=pbb(end-1,ii);
%             end
%             

            
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
%    end
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
