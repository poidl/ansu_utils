function regions=find_regions(vv)
% Connected component labeling: flood-fill approach
% see Steve Eddins' blog: http://blogs.mathworks.com/steve/2007/04/15/connected-component-labeling-part-4/
user_input;
[yi,xi]=size(vv);

% flag wet points (mixed layer excluded)
wet=~isnan(vv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative code

remove_south=[1:yi:xi*yi]-1; % southern neighbours of southern boundary points
remove_north=[yi:yi:xi*yi]+1; % northern neighbours of northern boundary points
%remove=[remove_south, remove_north];

L=zeros(size(wet)); % label matrix

iregion=1; % region label

while 1
    ii=find(wet==1,1,'first');
    if isempty(ii); break; end
    idx=ii; % linear indices of the pixels that were just labeled
    wet(idx)=0; % set the pixels to 0
    L(idx) = iregion; % label first point
    while ~isempty(idx); % find neighbours
        offsets = [-1; yi; 1; -yi]; % neighbor offsets for a four-connected neighborhood
        neighbors_tmp = bsxfun(@plus, idx, offsets'); % find all the nonzero 4-connected neighbors
        neighbors=setdiff(neighbors_tmp(:,1)',remove_south);
        neighbors=[neighbors, setdiff(neighbors_tmp(:,3)',remove_north)];
        neighbors=[neighbors, neighbors_tmp(:,2)', neighbors_tmp(:,4)'];
        if zonally_periodic;  
            neighbors(neighbors<1)=neighbors(neighbors<1)+xi*yi; 
            neighbors(neighbors>xi*yi)=neighbors(neighbors>xi*yi)-xi*yi;
        else
            neighbors=neighbors(neighbors>0 & neighbors<=xi*yi);
        end
        neighbors = unique(neighbors(:)); % remove duplicates
        idx = neighbors(wet(neighbors)); % keep only the neighbors that are nonzero
        if isempty(idx); break; end
        L(idx) = iregion;
        wet(idx)=0; % set the pixels that were just labeled to 0
    end
    iregion=iregion+1;
end

for ii=1:iregion-1;
    regions{ii}=find(L==ii);
end

% alternative code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative code

% cc=bwconncomp(wet,4);
% 
% if zonally_periodic; % in a zonally periodic domain merge regions which are cut apart by grid boundary
%     
%     % find regions which have points at western and eastern boundary
%     bdy_wet=false(1,length(cc.PixelIdxList));
%     ireg=1:length(cc.PixelIdxList);
%     for ii=ireg
%         if any(cc.PixelIdxList{ii}<=yi | cc.PixelIdxList{ii}>yi*(xi-1))
%             bdy_wet(ii)=true;
%         end
%     end
%     iw=ireg(bdy_wet);
%     
%     merged=false(1,length(cc.PixelIdxList));
%     
%     ii=1;
%     while ii<=length(iw)
%         for jj=1: length(iw)
%             if ii~=jj && ~(merged(iw(ii)) || merged(iw(jj)))
%                 pts1=cc.PixelIdxList{iw(ii)};
%                 pts2=cc.PixelIdxList{iw(jj)};
%                 % check if western border of region iw(ii) intersects
%                 % with eastern border of region iw(jj), or vice versa
%                 cond1=~isempty( intersect( pts1(pts1<=yi)+yi*(xi-1),pts2(pts2>yi*(xi-1)) ) );
%                 cond2=~isempty( intersect( pts2(pts2<=yi)+yi*(xi-1),pts1(pts1>yi*(xi-1)) ) );
%                 if cond1 | cond2
%                     cc.PixelIdxList{iw(ii)}=union( pts1, pts2 );
%                     merged(iw(jj))=true; % iw(jj) has been merged into iw(ii); delete later
%                     ii=ii-1;
%                     break
%                 end
%             end
%         end
%         ii=ii+1;
%     end
%     % delete merged regions
%     remove=ireg(merged);
%     for ii= remove(end:-1:1)
%         cc.PixelIdxList(ii)=[];
%     end
% end
% regions=cc.PixelIdxList;

% alternative code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


end 
