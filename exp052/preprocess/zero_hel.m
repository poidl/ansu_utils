% create ocean with zero helicity
clear all
close all
infile='data/input_data.mat';
load(infile);

ii=find(longs(1,1,:)>337,1,'first');

imagesc(lats(1,:,1),p(:,1,1),squeeze(ct(:,:,ii)));

d=~isnan(ct);
d=squeeze(sum(d,1));

if 0
    figure()
    imagesc(longs(1,1,:),lats(1,:,1),d)
    set(gca,'Ydir','normal')
    colorbar()
end


ct2=repmat(ct(:,:,ii),[1 1 size(longs,3)]);
s2=repmat(s(:,:,ii),[1 1 size(longs,3)]);

ct2(isnan(ct))=nan;
s2(isnan(ct))=nan;

if 0
    figure()
    imagesc(longs(1,1,:),lats(1,:,1),squeeze(ct2(20,:,:)))
    set(gca,'Ydir','normal')
    colorbar()
    figure()
    imagesc(longs(1,1,:),lats(1,:,1),squeeze(s2(20,:,:)))
    set(gca,'Ydir','normal')
    colorbar()
end

ct=ct2;
s=s2;

outfile=infile;
vars = {'s','ct','p','lats','longs'};
save(outfile,vars{:})
