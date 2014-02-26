function [drho,res]=solve_lsqr(regions, xx, yy)
user_input;
[yi,xi]=size(xx);
drho = nan(yi,xi);
load('data/stationindex.mat') % istation

for iregion=1:length(regions)
    
    region=regions{iregion};
    
    % set up east-west equations for weighted inversion
    reg=false(1,xi*yi)';
    reg(region)=true;
    
    en= reg & circshift(reg,-yi); %  find points between which a zonal gradient can be computed. en is true at a point if its eastward neighbor is in the region
    if ~zonally_periodic;  % remove equations for eastern boundary for zonally-nonperiodic domain
        en((xi-1)*yi+1:xi*yi)=false;
    end
    sreg=cumsum(reg); % sparse indices of points forming the region (points of non-region are indexed with dummy)
    sreg_en=circshift(sreg,-yi); % sparse indices of eastward neighbours
    
    j1_ew=sreg_en(en);  % j1 are j-indices for matrix coefficient 1
    j2_ew=sreg(en); % j2 are j-indices for matrix coefficient -1
    
    % set up north-south equations for weighted inversion
    nn= reg & circshift(reg,-1);
    nn(yi:yi:xi*yi)=false; % remove equations for northern boundary
    sreg_nn=circshift(sreg,-1);
    
    j1_ns=sreg_nn(nn);
    j2_ns=sreg(nn);
    
    % make the average of the potential zero
    % this should keep the surface from drifting away from the initial condition
    % we might change that to a different condition
    if ismember(istation,region)
        j1_condition=sreg(istation);
    else
        j1_condition=[1:sum(reg)];
    end
        
    j1=[j1_ew',j1_ns',j1_condition];
    j2=[j2_ew',j2_ns'];
    
    i2=1:(sum(en)+sum(nn)); % i-indices for matrix coeff. -1
    if ismember(istation,region)
        i1=[i2, (sum(en)+sum(nn)+1)]; % i-indices for matrix coeff. 1
    else
        i1=[i2, (sum(en)+sum(nn)+1)*ones(1,sum(reg))]; % i-indices for matrix coeff. 1
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new error measure
    inew=i1(end)+[1:sum(reg(:))];
    jnew=1:sum(reg(:));
    
    % build sparse matrices
    A=sparse([i1,i2,inew],[j1,j2,jnew],[ones(1,length(i1)), -ones(1,length(i2)), 0*ones(1,length(inew))]);    
    b=sparse( [xx(en); yy(nn); 0; zeros(sum(reg(:)),1) ]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    disp(['solving for region ',int2str(iregion)]);
    switch solver
        case 'iterative'
            %dbstop in lsqr at 321      
            [x,flag,relres] = lsqr(A,b,1e-15,5000);
            res=relres;
            if flag ~=0 
                disp('LSQR did not converge')
                if flag~=1
                    disp('MAXIT WAS NOT REACHED')
                end
            end
        case 'exact'
            x = (A'*A)\(A'*b);
    end
    
    x = full(x)';
    
    % put density changes calculated by the least-squares solver into
    % their appropriate position in the matrix

    drho(region) = x;
    
end
end
