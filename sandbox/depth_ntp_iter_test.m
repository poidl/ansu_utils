nz=5;
s1=linspace(35,36,nz)';
s2=linspace(36,37,nz)';
s=[s1 s2];

t1=linspace(11,10,nz)';
t2=linspace(12,13,nz)';
t=[t1 t2];

p=linspace(0,100,nz)';
p=[p p];

s0=[35.5 36.5];
t0=[10.5 12.5];
p0=[20 70];

[sns,tns,pns] = depth_ntp_iter(s0,t0,p0,s,t,p)





