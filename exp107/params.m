% TODO: the following should be set somewhere else
initial_surface_at_constant_pressure=false;
if initial_surface_at_constant_pressure;
    initial_pressure=1e3;
    nlevels=length(initial_pressure);
else
    glevels=[27.4];
    p_r=0; % reference pressure for initial surface
    nlevels=length(glevels);
end
