% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

runs=571:575;

fname=['../exp',num2str(runs(1)),'/data/iteration_history.mat'];
load(fname,'pns_hist');

nit=size(pns_hist,1);
[yi,xi]=size(squeeze(pns_hist(1,:,:)));

mask=true(yi,xi);

for ii=runs
    fname=['../exp',num2str(ii),'/data/iteration_history.mat'];
    load(fname,'pns_hist');
    pns=squeeze(pns_hist(nit,:,:));
    mask(isnan(pns))=false;
end


for ii=runs
    fname=['../exp',num2str(ii),'/data/mask.mat'];
    save(fname,'mask');
end



    
    
    


