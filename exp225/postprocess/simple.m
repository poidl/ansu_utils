clear all
load iteration_history
s1=sns_hist(1,:,1);
s2=sns_hist(2,:,1);
ct1=ctns_hist(1,:,1);
ct2=ctns_hist(2,:,1);
p1=pns_hist(1,:,1);
p2=pns_hist(2,:,1);



pmid1=0.5*(p1(1)+p1(2));
drhoy=gsw_rho(s1(2),ct1(2),pmid1)-gsw_rho(s1(1),ct1(1),pmid1);
gsw_rho(s2(2),ct2(2),pmid1)-gsw_rho(s2(1),ct2(1),pmid1)
