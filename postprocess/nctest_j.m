clear all
close all
f1='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp229/data/j1.nc';
f2='/home/nfs/z3439823/eclipse/workspace/ansu/j1.nc';
f3='/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp229/data/b.nc';
f4='/home/nfs/z3439823/eclipse/workspace/ansu/y.nc';
j1m=ncread(f1,'j1');
j1f=ncread(f2,'j1');
ym=ncread(f3,'b');
yf=ncread(f4,'y');

max(j1m)
max(j1f)

for i=1:length(j1m);
    indsm=find(j1m==j1m(i));
    if length(indsm)==2;
        ym_=ym(indsm);        
        indsf=find(j1f==j1m(i));
        yf_=yf(indsf);
        i
        
        if ~( all(ym_==yf_) || all(ym_([2 1])==yf_));
            error('hoit')
        end
    
    end
end