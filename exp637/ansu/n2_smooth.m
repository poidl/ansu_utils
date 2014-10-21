function [n2s,pmid]=n2_smooth(s,ct,p)

[nz,ny,nx]=size(s);
[n2,pmid]=gsw_Nsquared(s,ct,p);
n2=reshape(n2,[nz-1,ny,nx]);

n2s=repmat(nanmean(n2,3),[1,1,nx]);

if min(n2s(:)<=1e-8)
    error('problem')
end

% sloppy
n2s=cat(1,n2s(1,:,:),n2s);

end