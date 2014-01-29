% TODO: the following should be set somewhere else
initial_surface_at_constant_pressure=true;
%initial_surface_at_constant_pressure=true;
if initial_surface_at_constant_pressure;
    initial_pressure=1980.0;
    nlevels=length(initial_pressure);
else 
    glevels=[27.5];
    p_r=0; % reference pressure for initial surface
    nlevels=length(glevels);
end

