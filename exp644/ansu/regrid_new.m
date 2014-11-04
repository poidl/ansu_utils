function aout=regrid_new(a,dim,stag)
omega_user_input;
% stag
% - represents input grid, x represents output grid
%0:--- 1:- - -  (-1): - - - (2):- - - (-2): - - -  (99): none of the above
%  xxx    x x x      x x x       x x       x x x x


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(size(a))==3

    if dim==2
        a=permute(a,[1,3,2]);
    elseif dim==1
        a=permute(a,[2,3,1]);
    end

    if stag==-2
        aout=grd_append(a,3);
        aout=0.5*(aout+circshift(aout,[0 0 1]));
        in= isnan(aout) & ~isnan(a);
        aout(in)=a(in);
        if ~(dim==3 && zonally_periodic)
            aout(:,:,1)=a(:,:,1);
        end
    end


    if stag==2
        aout=0.5*(a+circshift(a,[0 0 -1]));
        if ~(dim==3 && zonally_periodic)
            aout(:,:,1)=a(:,:,1);
        end
    end

    if dim==2
        aout=permute(aout,[1,3,2]);
    elseif dim==1
        aout=permute(aout,[3,1,2]);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif length(size(a))==2

    if dim==1
        a=permute(a,[2,1]);
    end

    if stag==-2
        aout=grd_append(a,2);
        aout=0.5*(aout+circshift(aout,[0 1]));
        in= isnan(aout) & ~isnan(a);
        aout(in)=a(in);
        if ~(dim==2 && zonally_periodic)
            aout(:,1)=a(:,1);
        end
    end


    if stag==2
        aout=0.5*(a+circshift(a,[0 -1]));
        if ~(dim==2 && zonally_periodic)
            aout(:,1)=a(:,1);
        end
    end

    if dim==1
        aout=permute(aout,[2,1]);
    end

else
    error('error')
end
