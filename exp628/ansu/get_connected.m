function [sns,ctns,pns] = get_connected(sns,ctns,pns,istation)
    user_input;
    setnan=true(size(sns));
    regions=find_regions(sns);
    for iregion=1:length(regions)
        region=regions{iregion};
        if ismember(istation,region)
            setnan(region)=false;        
        end
    end 
    pns(setnan)=nan;
    sns(setnan)=nan;
    ctns(setnan)=nan;
    if no_land_mask 
        %load('data/no_land_mask.mat')
        error('this doesn''t work yet')
    end
end
