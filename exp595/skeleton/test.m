p1=ncread('data/nc/pns3d/pns3d1.nc','pns3d');
p2=ncread('data/nc/pns3d/pns3d2.nc','pns3d');

p1=permute(p1,[3 2 1]);
p2=permute(p2,[3 2 1]);

ps=cat(1,p1,p2);

checkit(p1)
checkit(p2)
checkit(ps) 

nk=size(p1,1);

ps=p1;
for ii=2:72
    p1=ncread(['data/nc/pns3d/pns3d',num2str(ii),'.nc'],'pns3d');
    p1=permute(p1,[3 2 1]);
    ps=cat(1,ps,p1);
    disp(ii);
    nk=[nk,size(p1,1)];
end
nk=cumsum(nk);
keyboard