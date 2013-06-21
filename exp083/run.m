if 1
    % add library paths
    restoredefaultpath
    addpath(genpath('../../../gsw_matlab_v3_02'))
    addpath(genpath('../external_scripts'))
    addpath(genpath('.'))

    % subsample data
    %subsample_gk_ak_gamma;

    % construct zero-helicity ocean from the subsampled data
    zero_hel_no_continents
    

    % run Paul's script
    omega_ak

    % perturb
    perturb;

    % run Paul's script with starting surface from above
    omega_ak_perturbed
end

%plot
 pressure_ns
 delta_pressure_ns
 phi_prime
 error_measure
%epsilon_vectors




