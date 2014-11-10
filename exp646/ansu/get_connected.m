function [sns,ctns,pns] = get_connected(sns,ctns,pns,istation)
    omega_user_input;
    setnan=true(size(sns));

    if no_land_mask 
        % dummy ex and ey
        ex=circshift(sns,[0 -1])-sns;
        ey=circshift(sns,[-1 0])-sns;
        if ~zonally_periodic
            ex(:,end)=nan;
        end
        ey(end,:)=nan;
        [ex,ey]=no_land_mask_disconnect(ex,ey);
        regions=find_regions_coupled_system(sns,ex,ey);
    else
        regions=find_regions(sns);
    end
    for iregion=1:length(regions)
        region=regions{iregion};
        if ismember(istation,region)
            setnan(region)=false;        
            pns(setnan)=nan;
            sns(setnan)=nan;
            ctns(setnan)=nan;
        end
    end 

end
