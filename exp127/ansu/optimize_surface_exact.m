function [sns_i,ctns_i,pns_i] = optimize_surface_exact(s,ct,p,sns,ctns,pns,e1t,e2t)

%           Optimize density surfaces to minimise the fictitious diapycnal diffusivity
%
%           Optimize density surface using an iterative method to minimise
%           fictitious diapycnal diffusivity according 'A new method of
%           forming approximately neutral surfaces, Klocker et al., Ocean
%           Science, 5, 155-172, 2009.'
%
%
% Input:    s           salinity
%           ct          conservative temperature
%           p           pressure
%           sns         salinity on initial density surface
%           ctns        conservative temperature on initial density
%                       surface
%           pns         pressure on initial density surface
%           e1t         grid length in meters in x-direction
%           e2t         grid length in meters in y-direction
%
% Output:   sns_i       salinity on optimized surface
%           ctns_i      conservative temperature on optimized surface
%           pns_i       pressure on optimized surface
%
% Calls:    mld.m, epsilon.m, var_on_surf.m
%
% Units:    salinity                    psu (IPSS-78)
%           conservative temperature    degrees C (IPS-90)
%           pressure                    dbar
%           gravitational acceleration  m/s^2
%           buoyancy frequency          s^-1
%           scale factors               m
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   Partially modified by P. Barker (2010-13)
%   Partially modified by S. Riha (2013)
%   Principal investigator: Trevor McDougall
%   type 'help analyze_surface' for more information
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

user_input;

% prepare data
%dbstop in mld at 50
cut_off_choice = mld(s,ct,p); % mixed-layer depth

stop_wetting=false;
% iterations of inversion
it=0; % start with it=0 and increment it after the initial surface is written to output
while it<=nit;
    
    % diagnose
    if save_iterations;
        if it==0; % dummy values
            phiprime=nan; ex=nan; ey=nan;
        end
        diagnose_and_write(it,sns,ctns,pns,ex,ey,phiprime);
    end
    
    if it==nit; % break out if done
        break
    end
    
    it=it+1; % start the next iteration
    disp(['Iteration nr.',int2str(it)]);
    
    
    % Locations where outcropping occurs may have changed. Add points to
    % surface if necessary.
%     if it<nit && ~stop_wetting ;
%         disp('Wetting')
%         if it==1;
%             [sns,ctns,pns,neighbours]=wetting(sns,ctns,pns,s,ct,p);
%         else
%             neighbours_old=neighbours;
%             [sns,ctns,pns,neighbours]=wetting(sns,ctns,pns,s,ct,p);
%             if all(neighbours==neighbours_old);
%                 stop_wetting=true;
%                 disp(['Wet points added are identical to previous iteration.'])
%             end
%         end
%     end
    
    % calculate slope errors/density gradient errors
    [ex,ey] = epsilon(sns,ctns,pns,e1t,e2t); 
    
    % disregard data above mixed layer depth
    ex(pns<= cut_off_choice)=nan;
    ey(pns<= cut_off_choice)=nan;
    sns(pns<=cut_off_choice)=nan;
    ctns(pns<=cut_off_choice)=nan;
    pns(pns<=cut_off_choice)=nan;
 
    
    % find independent regions -> a least-squares problem is solved for
    % each of these regions
    regions=find_regions(pns);
    
    % solve for phiprime
    phiprime=solve_lsqr(regions, ex, ey, e1t, e2t);
    
    % find corrected surface
    myfac=0;
    if myfac~=0;
        [sns_tmp, ctns_tmp, pns_tmp] = dz_from_phiprime(sns, ctns, pns, s, ct, p, phiprime );
        drho=myfac*(-2e-5*(pns_tmp-pns).*(ctns_tmp-ctns));
        [sns, ctns, pns] = dz_from_phiprime(sns_tmp, ctns_tmp, pns_tmp, s, ct, p, drho*1e-3 );
    else
        [sns, ctns, pns] = dz_from_phiprime(sns, ctns, pns, s, ct, p, phiprime );
    end
    
end

sns_i=sns;
ctns_i=ctns;
pns_i=pns;

end


function [sns,ctns,pns,neighbours]=wetting(sns,ctns,pns,s,ct,p)

[yi,xi]=size(sns);

wet=~isnan(squeeze(s(1,:,:))) & isnan(sns); % wet points at ocean surface excluding ans
wets=~isnan(sns); % wet points on ans
nn=wet(:) & circshift(wets(:),-1); % wet points with north. neighbour on ans
sn=wet(:) & circshift(wets(:),1);
en=wet(:) & circshift(wets(:),-yi);
wn=wet(:) & circshift(wets(:),yi);

% if a point adjacent to ans boundary has multiple neighbours, just do one neutral
% calculation
% TODO: start in eastward direction? Trevor has preference, see notes.
wn=wn & ~en;
nn=nn & ~wn & ~en;
sn=sn & ~nn & ~wn & ~en;

neighbour=circshift(en,yi);
[sns(en),ctns(en),pns(en)] = depth_ntp_iter(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,en),ct(:,en),p(:,en)); 

neighbour=circshift(wn,-yi);
[sns(wn),ctns(wn),pns(wn)] = depth_ntp_iter(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,wn),ct(:,wn),p(:,wn)); 

neighbour=circshift(nn,1);
[sns(nn),ctns(nn),pns(nn)] = depth_ntp_iter(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,nn),ct(:,nn),p(:,nn)); 

neighbour=circshift(sn,-1);
[sns(sn),ctns(sn),pns(sn)] = depth_ntp_iter(sns(neighbour)',ctns(neighbour)',pns(neighbour)',s(:,sn),ct(:,sn),p(:,sn)); 

neighbours= en | wn | nn | sn;

s1=sum(~isnan(sns(en)));
s2=sum(~isnan(sns(wn)));
s3=sum(~isnan(sns(nn)));
s4=sum(~isnan(sns(sn)));
disp(['Number of points added: ',num2str(s1+s2+s3+s4)])
end




function phiprime=solve_lsqr(regions, xx, yy, e1t, e2t)
user_input;
[yi,xi]=size(xx);
phiprime = nan(yi,xi);

for iregion=1:length(regions)
    
    region=regions{iregion};
    
    % set up east-west equations for weighted inversion
    reg=false(1,xi*yi)';
    reg(region)=true;
    
    en= reg & circshift(reg,-yi); %  find points between which phi_x can be computed. en is true at a point if its eastward neighbor is in the region
    if ~zonally_periodic;  % remove equations for eastern boundary for zonally-nonperiodic domain
        en((xi-1)*yi:xi*yi)=false;
    end
    sreg=cumsum(reg); % sparse indices of points forming the region (points of non-region are indexed with dummy)
    sreg_en=circshift(sreg,-yi); % sparse indices of eastward neighbours
    
    j1_ew=sreg_en(en);  % j1 are j-indices for matrix coefficient 1
    j2_ew=sreg(en); % j2 are j-indices for matrix coefficient -1
    
    % set up north-south equations for weighted inversion
    nn= reg & circshift(reg,-1);
    nn(1:yi:(xi-1)*yi+1)=false; % remove equations for northern boundary
    sreg_nn=circshift(sreg,-1);
    
    j1_ns=sreg_nn(nn);
    j2_ns=sreg(nn);
    
    % make the average of Phi' zero
    % this should keep the surface from drifting away from the initial condition
    % we might change that to a different condition
    j1_condition=[1:sum(reg)];
    
    j1=[j1_ew',j1_ns',j1_condition];
    j2=[j2_ew',j2_ns'];
    
    i2=1:(sum(en)+sum(nn)); % i-indices for matrix coeff. -1
    i1=[i2, (sum(en)+sum(nn)+1)*ones(1,sum(reg))]; % i-indices for matrix coeff. 1
    
    % build sparse matrices
    A=sparse([i1,i2],[j1,j2],[ones(1,length(i1)),-ones(1,length(i2))]);
    b=sparse( [xx(en).*e1t(en); yy(nn).*e2t(nn); 0 ]);
    
    disp(['solving for region ',int2str(iregion)]);
    switch solver
        case 'iterative'
            [x,dummyflag] = lsqr(A,b,1e-7,50000);
        case 'exact'
            x = (A'*A)\(A'*b);
    end
    
    x = full(x)';
    
    % put density changes calculated by the least-squares solver into
    % their appropriate position in the matrix

    phiprime(region) = x;
    
end
end


function [sns_out,ctns_out,pns_out] = dz_from_phiprime(sns, ctns, pns, s, ct, p, phiprime_e );

[zi,yi,xi]=size(s);
phiprime_e = permute(phiprime_e, [3 1 2]);

r=1.0;
delta = 1e-9;

rho_surf=gsw_rho(sns(:),ctns(:),pns(:));
t2=rho_surf.*(1-r*phiprime_e(:)); % approximately t2.*exp(-r*phiprime_e)

inds=1:yi*xi;
fr=true(1,yi*xi);
pns_out = nan(yi,xi);
sns_out = nan(yi,xi);
ctns_out = nan(yi,xi);
pns=pns(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative code

% stack=zi;
% ii=bsxfun(@times,1:yi*xi,ones(stack,1));
% pns_stacked=pns(ii); % stack pressure of current surface vertically
% t2_stacked_full=t2(ii); % stack locally referenced density of current surface vertically
% 
% t1=gsw_rho(s(:,:),ct(:,:),pns_stacked(:,inds)); % 3-d density referenced to pressure of the current surface
% t2_3d=t2_stacked_full(:,inds); %
% tni=t1-t2_3d; % rho-(rho_s+rho'); find corrected surface by finding roots of this term
% 
% Itni_p = tni>0;
% Itni_n = tni<0;
% 
% zc = any(Itni_p,1) & any(Itni_n,1); % horizontal indices of locations where zero-crossing occurs
% [min_tni, lminr]=min(abs(tni));
% cond1=min_tni>delta;
% final=min_tni<=delta; % at these horizontal indices root has been found
% fr= zc & cond1; %  at these horizontal locations we have to use root finding
% 
% lminr=lminr+stack*[0:size(Itni_n,2)-1];
% lminr=lminr(final);
% 
% sns_out(inds(final))=s(lminr); % adjust surface where root has already been found
% ctns_out(inds(final)) =ct(lminr);
% pns_out(inds(final)) =p(lminr);
% 
% k=sum(Itni_n,1); % find indices of flattened 3d-array where we use root finding
% k=k+stack*[0:size(Itni_n,2)-1];
%
% options=optimset('TolX',1e-3);
% for ii=inds(fr);
%     su=s(k(ii)); sl=s(k(ii)+1); % linear interpolation of s and ct
%     ctu=ct(k(ii)); ctl=ct(k(ii)+1);
%     pu=p(k(ii)); pl=p(k(ii)+1);
%     
%     sp=@(p) su+(sl-su)*(p-pu)/(pl-pu);
%     ctp=@(p) ctu+(ctl-ctu)*(p-pu)/(pl-pu);
%     fac=1./sqrt(gsw_rho(su,ctu,pu)^2+gsw_rho(sl,ctl,pl)^2); % scale to avoid having to set tolerance ?
%     drho_normalized=@(p)  fac.*(gsw_rho(sp(p), ctp(p), pns(ii)) -t2(ii));
%     
%     %dmyrho_normalized=@(p,su,sl,ctu,ctl,pu,pl)  1.*(gsw_rho( su+(sl-su)*(p-pu)/(pl-pu), ctu+(ctl-ctu)*(p-pu)/(pl-pu), pns(ii))-t2(ii));
%     %drho_normalized=@(p) dmyrho_normalized(p,su,sl,ctu,ctl,pu,pl);
%     proot=fzero(drho_normalized, [pu,pl],options);
%     
%     sns_out(ii)= su+(sl-su)*(proot-pu)/(pl-pu);
%     ctns_out(ii)=ctu+(ctl-ctu)*(proot-pu)/(pl-pu);
%     pns_out(ii)=proot;
% end

% alternative code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative code
    
refine_ints=100;

cnt=0;
while 1
    cnt=cnt+1;
    
    if cnt==1 % in first iteration pns_l is stacked vertically zi times, after that it is stacked refine_ints times
        stack=zi;
    elseif cnt==2
        stack=refine_ints+1;
    end
    if cnt==1 | cnt==2
        ii=bsxfun(@times,1:yi*xi,ones(stack,1));
        pns_stacked=pns(ii); % stack pressure of current surface vertically
        t2_stacked=t2(ii); % stack locally referenced density of current surface vertically
    end
    
    pns_stacked=pns_stacked(:,fr);
    t2_stacked=t2_stacked(:,fr);
    
    t1=gsw_rho(s(:,:),ct(:,:),pns_stacked); % 3-d density referenced to pressure of the current surface

    F=t1-t2_stacked; % rho-(rho_s+rho'); find corrected surface by finding roots of this term
    
    %dbstop in root_core at 11
    [final,fr,k_zc]=root_core(F,delta,stack);
    
    k_zc_3d=k_zc+stack*[0:size(F,2)-1]; % indices of flattened 3d-array where root has been found   
    
    sns_out(inds(final))=s(k_zc_3d(final)); % adjust surface where root has already been found
    ctns_out(inds(final)) =ct(k_zc_3d(final));
    pns_out(inds(final)) =p(k_zc_3d(final));
    inds=inds(fr); % points where surface has not been corrected
    
    if all(~fr) % break out of loop if all roots have been found
        break
    end
    
    k=k_zc_3d(fr);  % indices of flattened 3d-array where vertical resolution must be increased
    
    ds_ =  ( s(k+1) - s(k))/refine_ints; % increase resolution in the vertical
    dt_ = (ct(k+1) - ct(k))/refine_ints;
    dp_ =  (p(k+1) - p(k))/refine_ints;
    
    ds_ =bsxfun(@times, ds_, [0:refine_ints]');
    dt_ = bsxfun(@times, dt_, [0:refine_ints]');
    dp_ = bsxfun(@times, dp_, [0:refine_ints]');
    
    s =  bsxfun(@plus,s(k),ds_);
    ct =  bsxfun(@plus,ct(k),dt_);
    p =  bsxfun(@plus,p(k),dp_);

end

% alternative code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end




function regions=find_regions(vv)
% Connected component labeling: flood-fill approach
% see Steve Eddins' blog: http://blogs.mathworks.com/steve/2007/04/15/connected-component-labeling-part-4/
user_input;
[yi,xi]=size(vv);

% flag wet points (mixed layer excluded)
wet=~isnan(vv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative code

L=zeros(size(wet)); % label matrix

iregion=1; % region label

while 1
    ii=find(wet==1,1,'first');
    if isempty(ii); break; end
    idx=ii; % linear indices of the pixels that were just labeled
    wet(idx)=0; % set the pixels to 0
    L(idx) = iregion; % label first point
    while ~isempty(idx); % find neighbours
        offsets = [-1; yi; 1; -yi]; % neighbor offsets for a four-connected neighborhood
        neighbors = bsxfun(@plus, idx, offsets'); % find all the nonzero 4-connected neighbors
        if zonally_periodic;  
            neighbors(neighbors<1)=neighbors(neighbors<1)+xi*yi; 
            neighbors(neighbors>xi*yi)=neighbors(neighbors>xi*yi)-xi*yi;
        else
            neighbors=neighbors(neighbors>0 & neighbors<=xi*yi);
        end
        neighbors = unique(neighbors(:)); % remove duplicates
        idx = neighbors(wet(neighbors)); % keep only the neighbors that are nonzero
        if isempty(idx); break; end
        L(idx) = iregion;
        wet(idx)=0; % set the pixels that were just labeled to 0
    end
    iregion=iregion+1;
end

for ii=1:iregion-1;
    regions{ii}=find(L==ii);
end

% alternative code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative code

% cc=bwconncomp(wet,4);
% 
% if zonally_periodic; % in a zonally periodic domain merge regions which are cut apart by grid boundary
%     
%     % find regions which have points at western and eastern boundary
%     bdy_wet=false(1,length(cc.PixelIdxList));
%     ireg=1:length(cc.PixelIdxList);
%     for ii=ireg
%         if any(cc.PixelIdxList{ii}<=yi | cc.PixelIdxList{ii}>yi*(xi-1))
%             bdy_wet(ii)=true;
%         end
%     end
%     iw=ireg(bdy_wet);
%     
%     merged=false(1,length(cc.PixelIdxList));
%     
%     ii=1;
%     while ii<=length(iw)
%         for jj=1: length(iw)
%             if ii~=jj && ~(merged(iw(ii)) || merged(iw(jj)))
%                 pts1=cc.PixelIdxList{iw(ii)};
%                 pts2=cc.PixelIdxList{iw(jj)};
%                 % check if western border of region iw(ii) intersects
%                 % with eastern border of region iw(jj), or vice versa
%                 cond1=~isempty( intersect( pts1(pts1<=yi)+yi*(xi-1),pts2(pts2>yi*(xi-1)) ) );
%                 cond2=~isempty( intersect( pts2(pts2<=yi)+yi*(xi-1),pts1(pts1>yi*(xi-1)) ) );
%                 if cond1 | cond2
%                     cc.PixelIdxList{iw(ii)}=union( pts1, pts2 );
%                     merged(iw(jj))=true; % iw(jj) has been merged into iw(ii); delete later
%                     ii=ii-1;
%                     break
%                 end
%             end
%         end
%         ii=ii+1;
%     end
%     % delete merged regions
%     remove=ireg(merged);
%     for ii= remove(end:-1:1)
%         cc.PixelIdxList(ii)=[];
%     end
% end
% regions=cc.PixelIdxList;

% alternative code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


end 


function diagnose_and_write(it,sns,ctns,pns,ex,ey,phiprime_e)
user_input; % read nit, etc.

if it==0 % initialize
    [yi,xi]=size(sns);
    
    sns_hist = nan(nit+1,yi,xi); % store variables on initial surface and on nit improvements (=> nit+1)
    ctns_hist = nan(nit+1,yi,xi);
    pns_hist = nan(nit+1,yi,xi);
    
    slope_square = nan(nit,1);
    eps_rms_hist=nan(nit,1);
    phiprime_rms_hist=nan(nit,1);
    phiprime_e_hist = nan(nit,yi,xi);
    
    vars = {'sns_hist','ctns_hist','pns_hist','eps_rms_hist','phiprime_rms_hist','phiprime_e_hist'};
    save(history_file, vars{:},'-v7.3');
end

iteration_history = matfile(history_file,'Writable',true);
iteration_history.sns_hist(it+1,:,:) = permute(sns,[3 1 2]);
iteration_history.ctns_hist(it+1,:,:) = permute(ctns,[3 1 2]);
iteration_history.pns_hist(it+1,:,:) = permute(pns,[3 1 2]);

if it>0
    iteration_history.phiprime_e_hist(it,:,:) = permute(phiprime_e,[3,1,2]);
    
    s1=ex(~isnan(ex));  
    s2=ey(~isnan(ey));
    square=[s1(:) ; s2(:)].^2;
    slope_square(it,1) = sum(square);
    no_pts =length(square);
    iteration_history.eps_rms_hist(it,1) = sqrt(slope_square(it,1)/no_pts);
    
    iteration_history.phiprime_rms_hist(it,1)= sqrt( nanmean(phiprime_e(:).^2) );
    
end

end




