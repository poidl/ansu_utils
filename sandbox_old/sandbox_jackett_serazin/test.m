sa=[36. 36.]';
ct=[2. 1.]';
p=[0. 100.]';

out=nan(3,3);

tic
[out(1,1),out(1,2),out(1,3)]=depth_ntp(36., 1.5, 30., sa,ct,p);
display(['fzero   ',num2str(toc),' seconds']);
tic
[out(2,1),out(2,2),out(2,3)]=depth_ntp_serazin(36., 1.5, 30., sa,ct,p);
display(['serazin ',num2str(toc),' seconds']);
tic
[out(3,1),out(3,2),out(3,3)]=depth_ntp_jackett(36., 1.5, 30., sa,ct,p);
display(['jackett ',num2str(toc),' seconds']);

for ii=1:size(out,2);
    disp([num2str(out(1,ii)), '  ',num2str(out(2,ii)),' ',num2str(out(3,ii)) ])
end
    