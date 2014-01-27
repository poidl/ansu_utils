clear all;

fname='data/iteration_history.mat';
varname= 'sns_hist';
load(fname, varname);
varname= 'ctns_hist';
load(fname, varname);
varname= 'pns_hist';
load(fname, varname);

it=1;
p=1e3;
r_north=gsw_rho(sns_hist(it,end,1),ctns_hist(it,end,1),p);
r_south=gsw_rho(sns_hist(it,1,1),ctns_hist(it,1,1),p);


format long
dr=(r_north-r_south)

varname= 'drho_hist';
load(fname, varname);

dr1=drho_hist(it,1,1);
dr2=drho_hist(it,end,1);
(dr2+dr1)/dr2

