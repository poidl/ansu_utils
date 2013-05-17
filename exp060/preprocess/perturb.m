function perturb()


load('data/input_data.mat')

% location of perturbation
la=lats(1,:,1); lo=longs(1,1,:);
ilat=find(la>7,1,'first');
ilon=find(lo>200,1,'first');

ct=perturb_(ct,ilat,ilon);
s=perturb_(s,ilat,ilon);

vars = {'s','ct','p','lats','longs'};
save('data/input_data_perturbed.mat',vars{:})
end

function v=perturb_(v,ilat,ilon)

v(2:end,ilat,ilon)=v(1:end-1,ilat,ilon);
v(1,ilat,ilon)=v(2,ilat,ilon);

end