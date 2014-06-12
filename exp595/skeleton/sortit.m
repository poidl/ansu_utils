function [ilat,ilon,s,ct,p]=sortit(ilat,ilon,s,ct,p)

cnt=0;
while checkit(p)
    if cnt==5000
        disp('sorting takes long...')
        keyboard
    end
    [ilat,ilon,s,ct,p]=fixit(ilat,ilon,s,ct,p);
    cnt=cnt+1;   
end
disp(['finished sorting after ',num2str(cnt),' swaps'])

end