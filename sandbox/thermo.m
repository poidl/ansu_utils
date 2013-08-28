close all
clear all
% add library paths
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))

s=36.*ones(100,100);
ct=bsxfun(@times,linspace(-2,30,100),ones(100,1));
p=bsxfun(@times,ones(1,100),linspace(0,5000,100)');
alpha=gsw_alpha(s,ct,p);

contourf(p(:,1),ct(1,:), alpha)
xlabel('pressure [db]')
ylabel('ct [deg C]')
colorbar