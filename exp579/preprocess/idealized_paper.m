close all
clear all

xi=90;
yi=43;
zi=101;

ts=7; % surface temperature south
tn=5; % surface temperature north
dt=5; % surface temp minus bottom temp
dp=2e3;

y=linspace(0,1,yi);

tsn= ts+(tn-ts)*y;
tsn=repmat(tsn,[xi 1]);
tsn=repmat(permute(tsn,[3,2,1]),[zi 1 1]);

dt=-dt*linspace(0,1,zi);
dt=repmat(dt,[yi 1]);
dt=repmat(permute(dt,[2,1,3]),[1 1 xi]);

dp=dp*linspace(0,1,zi);
dp=repmat(dp,[yi 1]);
p=repmat(permute(dp,[2,1,3]),[1 1 xi]);

ct=tsn+dt;
s=0*ct+35;


contourf(squeeze(ct(:,:,1)))
colorbar()
figure()
rpot=gsw_rho(s,ct,0*ct);
contourf(squeeze(rpot(:,:,1)))
colorbar()

lats=[-84:4:84];
lats=repmat(lats,[xi 1]);
lats=repmat(permute(lats,[3,2,1]),[zi 1 1]);

longs=[0:4:356];
longs=repmat(longs,[yi 1]);
longs=repmat(permute(longs,[3,1,2]),[zi 1 1]);


b=1./(gsw_rho(s,ct,p).*gsw_alpha(s,ct,p));

vars = {'s','ct','p','lats','longs'};
%vars = {'s','ct','p'};
save('data/input_data.mat',vars{:})

save('data/b.mat','b')




