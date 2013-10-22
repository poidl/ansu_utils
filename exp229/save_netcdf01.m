function save_netcdf01(va,vname,fname);

nx=length(va);

delete(fname)

for ii=1:1;
    nccreate(fname,vname,...
              'Dimensions',{'pts' nx});
end
 
for ii=1:1;
    ncwrite(fname,vname, va);
end
