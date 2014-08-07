function save_netcdf_input_single_surf(sns,ctns,pns)

[ny,nx]=size(sns);

fname='data/os_input_single_surf.nc';


delete data/os_input_single_surf.nc


vname={'sns','ctns','pns'};

for ii=1:3;
    nccreate(fname,vname{ii},...
              'Dimensions',{'x' nx 'y' ny});
end

value={sns,ctns,pns};
 
for ii=1:3;
    ncwrite(fname,vname{ii}, permute(value{ii},[2 1]));
end

end

