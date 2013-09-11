restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))


rho_theta=1e3*gsw_alpha(36,10,1000)
rho_s=1e3*gsw_beta(36,10,1000)

D=4e3;
theta_z=10/D
s_z=5/D

kappa=1/gsw_sound_speed(35,10,1000)^2

left=rho_s*s_z+rho_theta*theta_z
right=kappa