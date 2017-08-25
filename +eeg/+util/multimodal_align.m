function [nirsout,eegout]=multimodalfusion(nirsin,eegin);
%This function aligns (with time-dillation) EEG and NIRS data based on
%common stimulus events

for idx=1:length(nirsin)
    
    x=lsqnonlin(@timedil,[0 0]);
    nirsout(idx)=nirsin(idx);    
    eegout(idx)=eegin(idx);
    eegout(idx).time=eegout(idx).time+x(1)+x(2)*eegout(idx).time;
    
end

    function cost=timedil(b)
        
        [X1, names1] = nirs.design.createDesignMatrix(nirsin(idx).stimulus,nirsin(idx).time);
        [X2, names2] = nirs.design.createDesignMatrix(eegin(idx).stimulus,nirsin(idx).time+b(1)+b(2)*nirsin(idx).time);
        X1=X1/max(X1(:));
        X2=X2/max(X2(:));
        
        lstBad1=~ismember(names1,names2);
        lstBad2=~ismember(names1,names2);
        X1(:,lstBad1)=[];
        X2(:,lstBad2)=[];
        cost=min(min(corrcoef(sum(X1,2),sum(X2,2))));
    end

[X1, names1] = nirs.design.createDesignMatrix(nirsout(1).stimulus,nirsin(1).time);
[X2, names2] = nirs.design.createDesignMatrix(eegout(1).stimulus,nirsin(1).time);
X1=X1/max(X1(:));
X2=X2/max(X2(:));
end