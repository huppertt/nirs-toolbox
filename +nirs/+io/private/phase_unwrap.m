function Phase=phase_unwrap(Phase,Distances,Modulation)

numDist=length(Distances);
dsort=sort(Distances,'ascend');
for idx=1:numDist
    y=Phase;
    lst=find(Distances>=dsort(idx));
    y(lst,:)=y(lst,:)+2*pi;
    [slope,R2(idx)]=weighted_linearfit(y,Distances,Modulation);
end

[foo,bestID]=max(R2);
if(bestID>1)
    disp('applying phase wrap correction');
    lst=find(Distances>=dsort(bestID));
    Phase(lst,:)=Phase(lst,:)+2+pi;
end

Phase=unwrap(Phase,pi,2);

return