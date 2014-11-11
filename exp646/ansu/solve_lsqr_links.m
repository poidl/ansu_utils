function [derr,res]=solve_lsqr_links(regions, xx, yy)
% does not assume that all neighbours have an equation linking them (contrary to solve_lsqr)
omega_user_input;

if clamp_on_point==true;
    load([datapath,'ibb.mat'])
end
[ny,nx]=size(xx);
derr = nan(ny,nx);

for iregion=1:length(regions)
    
    en=~isnan(xx(:));
    nn=~isnan(yy(:));
%     if ~zonally_periodic;  % remove equations for eastern boundary for zonally-nonperiodic domain
%         en((nx-1)*ny+1:nx*ny)=false;
%     end
%     nn(ny:ny:nx*ny)=false;
    
    region=regions{iregion};
    
    % set up east-west equations for weighted inversion
    reg=false(1,nx*ny)';
    reg(region)=true;
    
    en= en & reg;
    nn= nn & reg;
    
    sreg=cumsum(reg); % sparse indices of points forming the region (points of non-region are indexed with dummy)
    sreg_en=circshift(sreg,-ny); % sparse indices of eastward neighbours
    
    j1_ew=sreg_en(en);  % j1 are j-indices for matrix coefficient 1
    j2_ew=sreg(en); % j2 are j-indices for matrix coefficient -1
    
    sreg_nn=circshift(sreg,-1);
    
    j1_ns=sreg_nn(nn);
    j2_ns=sreg(nn);
    
    % clamp at point
    if ismember(ibb,region)
        j1_condition=sreg(ibb);
    else
        j1_condition=[1:sum(reg)];
    end
    
    j1=[j1_ew',j1_ns',j1_condition];
    j2=[j2_ew',j2_ns'];
    
    i2=1:(sum(en)+sum(nn)); % i-indices for matrix coeff. -1
    if ismember(ibb,region)
        i1=[i2, (sum(en)+sum(nn)+1)]; % i-indices for matrix coeff. 1
    else
        i1=[i2, (sum(en)+sum(nn)+1)*ones(1,sum(reg))]; % i-indices for matrix coeff. 1
    end
    
    % build sparse matrices
    A=sparse([i1,i2],[j1,j2],[ones(1,length(i1)),-ones(1,length(i2))]);
    b=sparse( [xx(en); yy(nn); 0 ]);
    
    %disp(['solving for region ',int2str(iregion)]);
    switch solver
        case 'iterative'
            tic
            [x,flag,relres,iter,resvec,lsvec] = lsqr(A,b,1e-15,10000);
            %display(['LSQR() took ',num2str(toc),' seconds for ',num2str(length(lsvec)),' iterations']);
            res=relres;
%             if flag ~=0 
%                 disp('LSQR did not converge')
%                 if flag~=1
%                     disp('MAXIT WAS NOT REACHED')
%                 end
%             end
        case 'exact'
            x = (A'*A)\(A'*b);
    end
    
    x = full(x)';
    
    % put density changes calculated by the least-squares solver into
    % their appropriate position in the matrix

    derr(region) = x;
    
end
end
