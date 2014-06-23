function flag=checkit(p1)
[nz,ny,nx]=size(p1);
di=diff(p1,1,1);
inds=di(:)<0;

flag=false;
if any(inds)
    flag=true;
    inds=di<0;
    for kk=1:size(inds,1)
        if any(inds(kk,:))
            disp(['inversion below level ',num2str(kk)])
%             ih=find(inds(kk,:));
%             for ll=ih
%                 [jj,ii]=ind2sub([ny,nx],ll);
%                 disp(['below level ',num2str(kk),': j=',num2str(jj),', i=',num2str(ii) ]);
%             end
        end
    end
end
end