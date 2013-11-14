clear all
close all

params; % load some parameters

load('data/input_data.mat')
for ii=1:size(s,3)
    for jj=1:size(s,2)
        kk=find(isnan(s(:,jj,ii)),1,'first');
        s(kk:end,jj,ii)=nan;
        ct(kk:end,jj,ii)=nan;
%        p(kk:end,jj,ii)=nan;
    end
end

lats=lats(1,:,1); longs=squeeze(longs(1,1,:))';
sa=s; clear s; % the _subs_ data saves sa in variable 's'

[longs,lats] = meshgrid(longs,lats);

[zi,yi,xi] = size(sa);

sns_i = nan(length(nlevels),yi,xi);
ctns_i = nan(length(nlevels),yi,xi);
pns_i = nan(length(nlevels),yi,xi);

Iak = 1;

if initial_surface_at_constant_pressure;
    ip=find(p(:,1,1)>=initial_pressure,1,'first');
    ptarget=p(ip,1,1)
    sns=squeeze(sa(ip,:,:));
    ctns=squeeze(ct(ip,:,:));
    pns=squeeze(p(ip,:,:));
    pns(isnan(sns))=nan;
else
    glevels_Iak=1e3+glevels(Iak);
    [zi,yi,xi]=size(sa);
    
    s_2=sa(:,:);
    ct_2=ct(:,:);
    p_2=p(:,:);
    inds=1:yi*xi;
    fr=true(1,yi*xi);
    pns_out = nan(yi,xi);
    sns_out = nan(yi,xi);
    ctns_out = nan(yi,xi);
    
    refine_ints=100;
    cnt=0;
    
    while 1
        cnt=cnt+1;
        
        if cnt==1 % in first iteration pns_l is stacked vertically zi times, after that it is stacked refine_ints times
            stack=zi;
        elseif cnt==2
            stack=refine_ints+1;
        end
        if cnt==1 | cnt==2
            ii=bsxfun(@times,1:yi*xi,ones(stack,1));
        end
        
        F=gsw_rho(s_2,ct_2,p_r*ones(size(s_2)))-glevels_Iak;
        [s_2,ct_2,p_2,sns_out,ctns_out,pns_out, inds, fr]=root_core(F,inds,refine_ints,s_2,ct_2,p_2,sns_out,ctns_out,pns_out);
        
        if all(~fr) % break out of loop if all roots have been found
            break
        end
        
    end
    sns=sns_out;
    ctns=ctns_out;
    pns=pns_out;
end


display('optimizing density surface');
tic
%dbstop in optimize_surface_exact at 624
%dbstop in optimize_surface_exact at 225
%dbstop in  optimize_surface_exact at 232 if ii==1220
%dbstop in optimize_surface_exact at 177
%dbstop in optimize_surface_exact at 104
%dbstop in optimize_surface_exact at 130


[sns_i(Iak,:,:),ctns_i(Iak,:,:),pns_i(Iak,:,:)] = optimize_surface_exact(sa,ct,p,sns,ctns,pns);
display(['optimizing density surface took ',num2str(toc),' seconds']);

        
        
  
