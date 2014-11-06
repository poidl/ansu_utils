function [sns,ctns,pns,nneighbours,iw]=wetting(sns,ctns,pns,s,ct,p)
% This function calls the actual wetting routine for each region
% separately, to avoid problems arising from two regions being separated by
% only one single wet point. In that case there could be wetting from only
% one of both directions (possibly the wrong one).
[yi,xi]=size(sns);

regions=find_regions(sns);
iw=false(xi*yi,1);
for ireg=1:length(regions)
    reg=regions{ireg};
    
    sns_r=nan*ones(yi,xi);
    ctns_r=nan*ones(yi,xi);
    pns_r=nan*ones(yi,xi);
    
    sns_r(reg)=sns(reg);
    ctns_r(reg)=ctns(reg);
    pns_r(reg)=pns(reg);
    

    [sns_r,ctns_r,pns_r,~,iw_r]=wetting_region(sns_r,ctns_r,pns_r,s,ct,p);

    sns(iw_r)=sns_r(iw_r);
    ctns(iw_r)=ctns_r(iw_r);
    pns(iw_r)=pns_r(iw_r);
    iw=iw|iw_r;
    
end

nneighbours=sum(iw);

end



function [sns,ctns,pns,nneighbours,iwetted]=wetting_region(sns,ctns,pns,s,ct,p)

    omega_user_input;

    [yi,xi]=size(sns);

    wet=~isnan(squeeze(s(1,:,:))) & isnan(sns); % wet points at ocean surface excluding ans
    wets=~isnan(sns); % wet points on ans
    nn=wet(:) & circshift(wets(:),-1); % wet points with north. neighbour on ans
    sn=wet(:) & circshift(wets(:),1);
    en=wet(:) & circshift(wets(:),-yi);
    wn=wet(:) & circshift(wets(:),yi);

    nn(yi:yi:yi*xi)=false;
    sn(1:yi:(xi-1)*yi+1)=false;
    if ~zonally_periodic;
        en((xi-1)*yi+1:xi*yi)=false;
        wn(1:yi)=false;
    end

    % if a point adjacent to ans boundary has multiple neighbours, just do one neutral
    % calculation
    % TODO: start in eastward direction? Trevor has preference, see notes.
    wn=wn & ~en;
    nn=nn & ~wn & ~en;
    sn=sn & ~nn & ~wn & ~en;

    inds=[1:xi*yi]';
    inds_neighbour=circshift(inds,-yi);
    neighbour=inds_neighbour(en);
    [sns(en),ctns(en),pns(en)] = depth_ntp_simple(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,en),ct(:,en),p(:,en)); 


    inds_neighbour=circshift(inds,yi);
    neighbour=inds_neighbour(wn);
    [sns(wn),ctns(wn),pns(wn)] = depth_ntp_simple(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,wn),ct(:,wn),p(:,wn)); 

    inds_neighbour=circshift(inds,-1);
    neighbour=inds_neighbour(nn);
    [sns(nn),ctns(nn),pns(nn)] = depth_ntp_simple(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,nn),ct(:,nn),p(:,nn)); 

    inds_neighbour=circshift(inds,1);
    neighbour=inds_neighbour(sn);
    [sns(sn),ctns(sn),pns(sn)] = depth_ntp_simple(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,sn),ct(:,sn),p(:,sn)); 

    s1=sum(~isnan(sns(en)));
    s2=sum(~isnan(sns(wn)));
    s3=sum(~isnan(sns(nn)));
    s4=sum(~isnan(sns(sn)));

    nneighbours=s1+s2+s3+s4;

    iwetted= ~isnan(sns(:)) & (en | wn | sn | nn);
    if sum(iwetted)~=nneighbours
        error('something is wrong')
    end
    %disp(['Number of points added: ',num2str(nneighbours)])

end



