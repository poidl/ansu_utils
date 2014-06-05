%[ocean, n] = gamma_ocean_and_n(sa,ct,p,lon,lat);
[nz,ny,nx]=size(sa);

land=isnan(squeeze(sa(1,:,:)));

sa_=36*ones(size(sa));
ct_=36*ones(size(sa));
[ocean, n] = gamma_ocean_and_n(sa_,ct_,p,lon,lat);
%keyboard
i1=ocean(:)==1 | ocean(:)==8; % North and Equatorial Pacific (west of Central America where land may be missing)
i2=~i1;



ct1=linspace(1,0,101)';
ct2=ct1+0.2;

ct1=repmat(ct1,[1 sum(i1)]);
ct2=repmat(ct2,[1 sum(i2)]);

ct(:,i1)=ct1;
ct(:,i2)=ct2;
sa(:,i1)=36;
sa(:,i2)=36.01;

land(25,70)=false; % remove central america
land(26,70)=false;
land(26,69)=false;
land(26,68)=false;

land(22,78)=true; % remove a (single-point) region so there's only one left

inan=false(ny,nx);
inan(:,1:35)=true;
inan(1:18,:)=true;
inan(35:end,:)=true;
inan(land(:))=true;
for kk=1:size(sa,1)
    sa(kk,inan(:))=nan;
    ct(kk,inan(:))=nan;
end

sa=reshape(sa,[nz,ny,nx]);
ct=reshape(ct,[nz,ny,nx]);
%keyboard