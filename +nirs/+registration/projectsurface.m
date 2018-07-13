function pos = projectsurface(pos,surf)

com = mean(surf,1);
for idx=1:size(pos,1)
    vec = pos(idx,:)-com;
     c = [0:.1:2*norm(vec)];
    vec=vec/norm(vec);
    p=c'*vec+ones(length(c),1)*com;
    [k,d]=dsearchn(surf,p);
    [~,i]=min(d);
    pos(idx,:)=p(i,:);
end

end