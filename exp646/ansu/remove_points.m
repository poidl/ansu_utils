function [regions]=remove_points(lnbx,lnby,pns)
    goodx=~isnan(lnbx);
    goody=~isnan(lnby);
    
    % lnbx and lnby are staggered with respect to pns
    % tmp lives on the pns grid and is true if at least one (staggered) neighbour is not-nan.
    tmp=goodx | circshift(goodx,[0 1]);
    tmp= tmp | (goody | circshift(goody,[1 0]));
    
    % eqe is true if a tmp has an eastern neighbour...
    eqe=tmp & circshift(tmp,[0 -1]);
    eqn=tmp & circshift(tmp,[-1 0]);
    
    % ...there is an eastern neighbour, but not an eastward equation. Remove.
    remove_eq= eqe&~goodx;
    remove_eq=remove_eq | eqn&~goody;
    
    tmp(remove_eq)=false;
    
    tmp=double(tmp);
    tmp(tmp==0)=nan;
    
    disp(['      Lost ', num2str(sum(~isnan(pns(:)))-sum(~isnan(tmp(:)))), ...
              ' points during re-gridding'])
    if any(tmp(:))
        regions=find_regions(tmp);
    else
        regions={};
    end
end
