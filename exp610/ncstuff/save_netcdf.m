function save_netcdf(va,vname,fname);

ss=size(va);
if length(ss)==3  
    cs=3;
    [nz,ny,nx]=size(va);
elseif length(ss)==2
    cs=2;
    [ny,nx]=size(va);
end

delete(fname)

if cs==3
    nccreate(fname,vname,...
              'Dimensions',{'x' nx 'y' ny 'z' nz});
    ncwrite(fname,vname, permute(va,[3 2 1]));
elseif cs==2
    nccreate(fname,vname,...
              'Dimensions',{'x' nx 'y' ny});
    ncwrite(fname,vname, permute(va,[2 1]));  
end
