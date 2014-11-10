
    %disp('Appending') 
    
    if ~stop_wetting
        [sns,ctns,pns,nneighbours,iwetted]=wetting(sns,ctns,pns,s,ct,p);
    else
        nneighbours=0;
    end
    for ll=1:length(iwetted_old)
        iw_old=iwetted_old{ll};
        if size(iwetted)==size(iw_old)
            app_and_cut=all(iwetted==iw_old); % identical points keep getting appended and cut off
            if app_and_cut
                stop_wetting=true;
            end
        end
    end
    iwetted_old(1)=[]; % forget last
    iwetted_old{1,3}=iwetted; % insert most recent
    if nneighbours==0
        stop_wetting=true;
    end
    if it>5  
        if stop_wetting
            cnt_it_after_wetting=cnt_it_after_wetting+1;
        end
    end
    if cnt_it_after_wetting==nit_after_wetting
        breakout=true;
    end
    
    
    