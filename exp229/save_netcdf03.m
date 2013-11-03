function save_netcdf03(va,vname,fname);

[nz,ny,nx]=size(va);

delete(fname)

for ii=1:1;
    nccreate(fname,vname,...
              'Dimensions',{'x' nx 'y' ny 'z' nz});
end
 
for ii=1:1;
    ncwrite(fname,vname, permute(va,[3 2 1]));
end