function p2=fixit(p1)

[nz,ny,nx]=size(p1);
di=diff(p1,1,1);
inds=di(:)<0;

if any(inds)
    inds=di<0;
    for kk=1:size(inds,1)
        if any(inds(kk,:))
            tmp=p1(kk,:,:);
            p1(kk,:,:)=p1(kk+1,:,:);
            p1(kk+1,:,:)=tmp;
        end
    end
end
p2=p1;

end
