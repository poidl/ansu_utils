% TODO: the following should be set somewhere else
initial_surface_at_constant_pressure=false;
%initial_surface_at_constant_pressure=true;
if initial_surface_at_constant_pressure;
    initial_pressure=1e3;
    nlevels=length(initial_pressure);
else
    % idealized_01:   
    % sigma0: mean pressure,  28.2: 887,  28.3: 1045,  28.25: 965,  28.275: 1004,   28.2719464:  1000.0005
    glevels=[28.2719464];
    p_r=0; % reference pressure for initial surface
    nlevels=length(glevels);
end
