function [ilat,ilon,s,ct,p]=fixit(ilat,ilon,s,ct,p)

[nz,ny,nx]=size(p);
di=diff(p,1,1);
inds=di(:)<0;

if any(inds)
    done=false;
    inds=di<0;
    for kk=1:size(inds,1)
        if any(inds(kk,:))
            tmp=ilat(kk);
            ilat(kk)=ilat(kk+1);
            ilat(kk+1)=tmp;
            
            tmp=ilon(kk);
            ilon(kk)=ilon(kk+1);
            ilon(kk+1)=tmp;
            
            tmp=s(kk,:,:);
            s(kk,:,:)=s(kk+1,:,:);
            s(kk+1,:,:)=tmp;

            tmp=ct(kk,:,:);
            ct(kk,:,:)=ct(kk+1,:,:);
            ct(kk+1,:,:)=tmp;
            
            tmp=p(kk,:,:);
            p(kk,:,:)=p(kk+1,:,:);
            p(kk+1,:,:)=tmp;
            
            disp(['... swapped levels ',num2str(kk), ' and ',num2str(kk+1)])
            done=true; % only swap once
        end
        if done
            break
        end
    end
end

end
