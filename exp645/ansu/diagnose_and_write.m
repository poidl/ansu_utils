function diagnose_and_write(it,sns,ctns,pns,er_delx,er_dely,derr,res,b,n2ns)
omega_user_input; % read nit, etc.

if it==0 % initialize
    [yi,xi]=size(sns);
    
    sns_hist = nan(1,yi,xi); 
    ctns_hist = nan(1,yi,xi);
    pns_hist = nan(1,yi,xi);
    
    %slope_square = nan(nit,1);
    %eps_rms_hist=nan(nit,1);
    
    err_del_rms_hist=nan(1,1); 
    derr_rms_hist=nan(1,1);
    err_grad_rms_hist=nan(1,1);
    
    derr_hist = nan(1,yi,xi);
    res_hist = nan(1,1);
    
    b_hist = nan(1,yi,xi);
    n2ns_hist = nan(1,yi,xi);
    
    vars = {'sns_hist','ctns_hist','pns_hist','err_del_rms_hist','derr_rms_hist','err_grad_rms_hist','derr_hist','res_hist','b_hist','n2ns_hist'};
    save(history_file, vars{:},'-v7.3');
end

iteration_history = matfile(history_file,'Writable',true);
iteration_history.sns_hist(it+1,:,:) = permute(sns,[3 1 2]);
iteration_history.ctns_hist(it+1,:,:) = permute(ctns,[3 1 2]);
iteration_history.pns_hist(it+1,:,:) = permute(pns,[3 1 2]);

if it>0
    iteration_history.derr_hist(it,:,:) = permute(derr,[3,1,2]);
    iteration_history.res_hist(it,1) = res;
    
    iteration_history.b_hist(it,:,:) = permute(b,[3,1,2]);
    iteration_history.n2ns_hist(it,:,:) = permute(n2ns,[3,1,2]);
    
    %s1=ex(~isnan(ex));  
    %s2=ey(~isnan(ey));
    %square=[s1(:) ; s2(:)].^2;
    %slope_square(it,1) = sum(square);
    %no_pts =length(square);
    %iteration_history.eps_rms_hist(it,1) = sqrt(slope_square(it,1)/no_pts);
    
    % not area weighted:
    iteration_history.derr_rms_hist(it,1)= sqrt( nanmean(derr(:).^2) );
    iteration_history.err_del_rms_hist(it,1)= sqrt( nanmean( [er_delx(:);er_dely(:)] .^2)); % staggerd grid
    load([datapath,'dxdy.mat'])
    er_gradx=er_delx./dx;
    er_grady=er_dely./dy;
    iteration_history.err_grad_rms_hist(it,1)= sqrt( nanmean( [er_gradx(:);er_grady(:)] .^2)); % staggerd grid
   
    [gradx_times_sqrtdA,grady_times_sqrtdA]=times_sqrtdA_on_delta(er_gradx,er_grady,dx,dy);

    myarea=dx.*dy;
    A=sum(myarea(~isnan(sns))); % total area    
    weighted=sqrt( (1./A)*nansum( [gradx_times_sqrtdA(:);grady_times_sqrtdA(:)].^2 ) );
    
    iteration_history.err_grad_rms_area_weighted_hist(it,1)= weighted;
    
end

end
