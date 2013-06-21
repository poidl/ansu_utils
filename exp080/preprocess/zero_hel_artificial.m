clear all
close all

% size of subs_4x4 is (101,43,90)
sz=[101,3,3];


lsdepth=linspace(0,1,sz(1))';
lsdepth=repmat( permute(lsdepth,[1,2,3]),[1,sz(2),sz(3)] );

lslat=linspace(0,1,sz(2))';
lslat=repmat(permute(lslat,[3,1,2]),[sz(1),1,sz(3)]);

% first two coefficients determine vertical gradient;  ct(end,end,1), ct(1,end,1),
% third is 0, fourth determines meridional gradient; (ct(1,1,1)-ct(1,end,1))
ct=0.*lsdepth+20.*(1-lsdepth)+0.*lslat+1.*(1-lslat);
s= 35.*lsdepth+36.*(1-lsdepth);


p=3000.*lsdepth;

lats=90.*lslat;

lslon=linspace(0,1,sz(3));
lslon=repmat(permute(lslon,[3,1,2]),[sz(1),sz(2),1]);
longs=360.*lslon;

vars = {'s','ct','p','lats','longs'};
save('data/input_data.mat',vars{:})

