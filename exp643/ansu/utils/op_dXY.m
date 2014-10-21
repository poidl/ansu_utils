function f=op_dXY(f,op,dXY,staggering)

% op: 'times' or 'over'
% dXY: 'dx','dy', 'dA' or 'mean_dxy'
% staggering='i', 'j' or 'none'
%   'i': f lives on the same j-grid as s,t,p, but is staggered in the i-dimension
%   'j': f lives on the same i-grid as s,t,p, but is staggered in the j-dimension
%   'none': no staggering

    load('data/dxdy.mat')

    if strcmp(dXY,'dx')
        dd=get_dx(dx,staggering);

    elseif strcmp(dXY,'dy')
        dd=get_dy(dy,staggering);

    elseif strcmp(dXY,'dA')
        dx=get_dx(dx,staggering);
        dy=get_dy(dy,staggering);
        dd=dx.*dy;

    elseif strcmp(dXY,'sqrtdA')
        dx=get_dx(dx,staggering);
        dy=get_dy(dy,staggering);
        dd=sqrt(dx.*dy);
        
    elseif strcmp(dXY,'mean_dxy')
        dx=get_dx(dx,staggering);
        dy=get_dy(dy,staggering);
        dd=0.5*(dx+dy);        

    else
        error('problem')
    end
    
    
    nxy=length(dx(:));

    if length(f(:))==nxy % f is probably a lateral surface
        if strcmp(op,'times')
            f(:)=f(:).*dd(:);
        elseif strcmp(op,'over')
            f(:)=f(:)./dd(:);
        else
            error('problem')
        end                
    elseif size(f(:,:),2)==nxy % f is probably a 3d field
        if strcmp(op,'times')    
            f(:,:)=bsxfun(@times,f(:,:),dd(:)');
        elseif strcmp(op,'over')        
            f(:,:)=bsxfun(@times,f(:,:),1./dd(:)'); 
        else
            error('problem')
        end
    else
        error('problem')
    end 
    
    

end



function dx=get_dx(dx,staggering)
    if strcmp(staggering,'i')
        % nothing to do
        
    elseif strcmp(staggering,'none')
        dx=0.5*(dx+circshift(dx,[0 1]));
        dx(:,1)=dx(:,2); % sloppy 
        
    elseif strcmp(staggering,'j')
        dx=0.5*(dx+circshift(dx,[0 1]));
        dx(:,1)=dx(:,2); % sloppy
        dx=0.5*(dx+circshift(dx,[-1 0]));
        dx(end,:)=dx(end-1,:); % sloppy  
        
    else
        error('problem')
    end
end


        
function dy=get_dy(dy,staggering)
    if strcmp(staggering,'j')
        % nothing to do
        
    elseif strcmp(staggering,'none')
        dy=0.5*(dy+circshift(dy,[1 0]));
        dy(1,:)=dy(2,:); % sloppy
        
    elseif strcmp(staggering,'i')
        dy=0.5*(dy+circshift(dy,[1 0]));
        dy(1,:)=dy(2,:); % sloppy
        dy=0.5*(dy+circshift(dy,[0 -1]));
        dy(:,end)=dy(:,end-1); % sloppy
        
    else
        error('problem')
    end
end        


% 
% function res=times(f,dXY)
% 
%     nxy=length(dx(:));
% 
%     if length(f(:))==nxy % f is probably a lateral surface
%         if strcmp(op,'times')
%             f(:)=f(:).*dA(:);
%         elseif strcmp(op,'over')
%             f(:)=f(:)./dA(:);
%         else
%             error('problem')
%         end                
%     elseif size(f(:,:),2)==nxy % f is probably a 3d field
%         if strcmp(op,'times')    
%             f(:,:)=bsxfun(@times,f(:,:),dA(:)');
%         elseif strcmp(op,'over')        
%             f(:,:)=bsxfun(@times,f(:,:),1./dA(:)'); 
%         else
%             error('problem')
%         end
%     else
%         error('problem')
%     end
%     
% end
