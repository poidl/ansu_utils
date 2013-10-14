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

phi=(r_north-r_south)/r_south

varname= 'phiprime_e_hist';
load(fname, varname);

phiprime_e_hist(it,1,1)-phiprime_e_hist(it,1,1)
