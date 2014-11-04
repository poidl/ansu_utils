function a=grd_append(a,dim)

if length(size(a))==3

    [nz,ny,nx]=size(a);

    if dim==2
        a=permute(a,[1,3,2]);
    elseif dim==1
        a=permute(a,[2,3,1]);
    end

    aw=circshift(a,[0 0 1]);
    in= isnan(a) & ~isnan(aw);
    a(in)=aw(in);


    if dim==2
        a=permute(a,[1,3,2]);
    elseif dim==1
        a=permute(a,[3,1,2]);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif length(size(a))==2

    [ny,nx]=size(a);

    if dim==1
        a=permute(a,[2,1]);
    end

    aw=circshift(a,[0 1]);
    in= isnan(a) & ~isnan(aw);
    a(in)=aw(in);

    if dim==1
        a=permute(a,[2,1]);
    end

else
    error('error')
end
