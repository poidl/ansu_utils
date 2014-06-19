close all

load('postprocess/bb_difference.mat')

va1=squeeze(ps(jsa,:,:));
va2=squeeze(ps(js,:,:));


sz=1.0*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 
va=va1;
h=imagesc(va);
set(h,'alphadata',~isnan(va))
set(gca,'YDir','normal')
colorbar()
hold on
plot(ilo,ila,'w*','markersize',30)
%plot(ilo_,ila_,'y*','markersize',30)
title('pressure of surface [dbar]')
print('-dpdf','-r200',['figures/bb_difference01.pdf'])

figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 
va=va2;
h=imagesc(va);
set(h,'alphadata',~isnan(va))
set(gca,'YDir','normal')
colorbar()
hold on
%plot(ilo,ila,'w*','markersize',30)
plot(ilo_,ila_,'y*','markersize',30)
title('pressure of surface [dbar]')
print('-dpdf','-r200',['figures/bb_difference02.pdf'])

figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 
di=va2-va1;
va=log10(abs(di));
h=imagesc(va);
set(h,'alphadata',~isnan(va))
set(gca,'YDir','normal')
colorbar()
hold on
plot(ilo,ila,'w*','markersize',30)
plot(ilo_,ila_,'y*','markersize',30)
title('|p1-p2| [log10(dbar)]')
print('-dpdf','-r200',['figures/bb_difference03.pdf'])


va1(ila,ilo)
va2(ila,ilo)
abs(va1(ila,ilo)-va2(ila,ilo))
va1(ila_,ilo_)
va2(ila_,ilo_)