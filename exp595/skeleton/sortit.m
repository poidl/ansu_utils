function p2=sortit(p1)

cnt=0;
while checkit(p1)
    if cnt==50
        disp('sorting takes long...')
        keyboard
    end
    p1=fixit(p1);
    cnt=cnt+1;   
end
disp(['finished sorting after ',num2str(cnt),' iterations'])
p2=p1;

end