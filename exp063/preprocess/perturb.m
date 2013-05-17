function perturb()


load('data/input_data.mat')

% location of perturbation
la=lats(1,:,1); lo=longs(1,1,:);
ilat=find(la>=45,1,'first')-1;
ilon=find(lo>=180,1,'first')-1;

ct=perturb_(ct,ilat,ilon);
s=perturb_(s,ilat,ilon);

vars = {'s','ct','p','lats','longs'};
save('data/input_data_perturbed.mat',vars{:})
end

function v=perturb_(v,ilat,ilon)

v(2:end,ilat,ilon)=v(1:end-1,ilat,ilon);
v(1,ilat,ilon)=v(2,ilat,ilon);

end