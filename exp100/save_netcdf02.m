function save_netcdf02(va,vname,fname);

[ny,nx]=size(va);

delete(fname)

for ii=1:1;
    nccreate(fname,vname,...
              'Dimensions',{'x' nx 'y' ny});
end
 
for ii=1:1;
    ncwrite(fname,vname, permute(va,[2 1]));
end
