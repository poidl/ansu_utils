% TODO: the following should be set somewhere else
initial_surface_at_constant_pressure=false;
%initial_surface_at_constant_pressure=true;
if initial_surface_at_constant_pressure;
    initial_pressure=1e3;
    nlevels=length(initial_pressure);
else % 27.5: 974  27.6: 1140  27.55: 1045   27.54 1029
    glevels=[27.537];
    p_r=0; % reference pressure for initial surface
    nlevels=length(glevels);
end
