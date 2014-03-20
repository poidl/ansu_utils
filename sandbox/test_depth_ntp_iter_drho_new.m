addpath(genpath('../exp583'))

p=0:20:100;
p=repmat(p,[2,1]);
p=repmat(p,[1 1 2]);
p=permute(p,[2 1 3]);

s=0*p+36;

ct=zeros(6,2,2);
n=6;
t1=5;
t2=6;
dct=10;
ct(:,1,1)=linspace(t1,t2,n);
ct(:,1,2)=linspace(t1,t2,n);
ct(:,2,1)=linspace(t1,t2,n);
ct(:,2,2)=linspace(t1+dct,t2+dct,n);

s0=36*ones(2,2);
ct0=0.5*(t1+t2)*ones(2,2);
p0=50*ones(2,2);

depth_ntp_iter_drho_new(s0,ct0,p0,s,ct,p,0*p);
