
% plot scalar as function of iterations, e.g. rms slope error
clear all;
close all;

runs=[444,550,571];
fname=['../exp',num2str(runs(1)),'/data/input_data.mat'];
load(fname);
[zi,yi,xi]=size(s);
bot1=nan*ones(yi,xi);
bot2=nan*ones(yi,xi);
bot3=nan*ones(yi,xi);

cnt=0;
for rr=runs
    cnt=cnt+1;
    fname=['../exp',num2str(rr),'/data/b.mat'];
    load(fname);
    fname=['../exp',num2str(rr),'/data/input_data.mat'];
    load(fname);
    [zi,yi,xi]=size(s);
    if cnt==1
        b1=b;
        lat1=squeeze(lats(1,:,1));
        lon1=squeeze(longs(1,1,:));
        p1=squeeze(p(:,1,1));
        for ii=1:xi
            for jj=1:yi
                ibot=find(~isnan(s(:,jj,ii)),1,'last');
                if ~isempty(ibot)
                    bot1(jj,ii)=p(ibot,jj,ii);
                end
            end
        end
    elseif cnt==2
        b2=b;
        lat2=squeeze(lats(1,:,1));
        lon2=squeeze(longs(1,1,:));
        p2=squeeze(p(:,1,1));
        for ii=1:xi
            for jj=1:yi
                ibot=find(~isnan(s(:,jj,ii)),1,'last');
                if ~isempty(ibot)
                    bot2(jj,ii)=p(ibot,jj,ii);
                end
            end
        end        
    elseif cnt==3
        b3=b;
        lat3=squeeze(lats(1,:,1));
        lon3=squeeze(longs(1,1,:));   
        p3=squeeze(p(:,1,1));        
        for ii=1:xi
            for jj=1:yi
                ibot=find(~isnan(s(:,jj,ii)),1,'last');
                if ~isempty(ibot)
                    bot3(jj,ii)=p(ibot,jj,ii);
                end
            end
        end             
    end
end


%sz=40*[1.2 1];
sz=40*[1.5 1];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 
set(gcf,'DefaultAxesFontSize', 18)
set(gcf,'DefaultTextFontSize',18)

lo=188;
[mi,l1]=min(abs(lon1-lo));
lo=330;
[mi,l2]=min(abs(lon1-lo));
disp(['lon1: ',num2str(lon1(l1))])
disp(['lon2: ',num2str(lon1(l2))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sub=subplot(3,2,1)
vp=squeeze(b1(:,:,l1));
h=imagesc(lat1,p1, vp)
set(h,'alphadata',~isnan(vp)) % white nans
hold on
plot(lat1,bot1(:,l1))
colorbar
p = get(sub, 'position');
p(1) = p(1)-0.05;
p(3) = 1.5*p(3);
set(sub, 'position', p);

sub=subplot(3,2,2)
vp=squeeze(b1(:,:,l2));
h=imagesc(lat1,p1, vp)
set(h,'alphadata',~isnan(vp)) % white nans
hold on
plot(lat1,bot1(:,l2))
colorbar
p = get(sub, 'position');
p(1) = p(1)-0.05;
p(3) = 1.5*p(3);
set(sub, 'position', p);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sub=subplot(3,2,3)
vp=squeeze(b2(:,:,l1));
h=imagesc(lat2,p2, vp)
set(h,'alphadata',~isnan(vp)) % white nans
hold on
plot(lat2,bot2(:,l1))
colorbar
p = get(sub, 'position');
p(1) = p(1)-0.05;
p(3) = 1.5*p(3);
set(sub, 'position', p);

sub=subplot(3,2,4)
vp=squeeze(b2(:,:,l2));
h=imagesc(lat2,p2, vp)
set(h,'alphadata',~isnan(vp)) % white nans
hold on
plot(lat2,bot2(:,l2))
colorbar
p = get(sub, 'position');
p(1) = p(1)-0.05;
p(3) = 1.5*p(3);
set(sub, 'position', p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sub=subplot(3,2,5)
vp=squeeze(b3(:,:,l1));
h=imagesc(lat3,p3, vp)
set(h,'alphadata',~isnan(vp)) % white nans
hold on
plot(lat3,bot3(:,l1))
colorbar
p = get(sub, 'position');
p(1) = p(1)-0.05;
p(3) = 1.5*p(3);
set(sub, 'position', p);

sub=subplot(3,2,6)
vp=squeeze(b3(:,:,l2));
h=imagesc(lat3,p3, vp)
set(h,'alphadata',~isnan(vp)) % white nans
hold on
plot(lat3,bot3(:,l2))
colorbar
p = get(sub, 'position');
p(1) = p(1)-0.05;
p(3) = 1.5*p(3);
set(sub, 'position', p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


print('-dpng','-r200',['figures/b'])
