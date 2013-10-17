clear all
load iteration_history
% s1=sns_hist(1,:,1);
% s2=sns_hist(2,:,1);
% ct1=ctns_hist(1,:,1);
% ct2=ctns_hist(2,:,1);
% p1=pns_hist(1,:,1);
% p2=pns_hist(2,:,1);
% 
% 
% 
% pmid1=0.5*(p1(1)+p1(2));
% drhoy=gsw_rho(s1(2),ct1(2),pmid1)-gsw_rho(s1(1),ct1(1),pmid1)
% gsw_rho(s2(2),ct2(2),pmid1)-gsw_rho(s2(1),ct2(1),pmid1)
% 
% pmid2=0.5*(p2(1)+p2(2));
% gsw_rho(s2(2),ct2(2),pmid2)-gsw_rho(s2(1),ct2(1),pmid2)


% first term in Eq. 18
sns=squeeze( sns_hist(:,2,1) );
ctns=squeeze( ctns_hist(:,2,1) );
pns=squeeze( ctns_hist(:,2,1) );

Delta_plus1=nan*ones(1,length(sns)-1);

for it=1:9;
    pmid1=0.5*(pns_hist(it,1,1)+pns_hist(it,2,1));
    %pmid2=0.5*(pns_hist(it+1,1,1)+pns_hist(it+1,2,1));

    rBpi=gsw_rho(sns(it),ctns(it),pmid1);
    rCpi=gsw_rho(sns(it+1),ctns(it+1),pmid1);

    rBpi_=gsw_rho(sns(it),ctns(it),pmid1+1);
    rCpi_=gsw_rho(sns(it+1),ctns(it+1),pmid1+1);

    Delta_plus1(it)=((rBpi_-rBpi)-(rCpi_-rCpi))*500;
end






