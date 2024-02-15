function dataout=generate_matching_dataset(n,sampledata,sampleStats)

if(nargin<3)
    j=nirs.modules.GLM;
    sampleStat=j.run(sampledata);
end

pMax = 4*sampledata(1).Fs;

for i=1:length(sampledata)
    mu{i}=mean(sampledata(i).data,1);
    [inn{i},f{i}]=nirs.math.innovations(sampledata(i).data-ones(size(sampledata(i).data,1),1)*mu{i},pMax);
end
inn=vertcat(inn{:});

for i=1:n
    dataout(i,1)=sampledata(randi(length(sampledata),1));
    lst=randperm(size(inn,1),length(dataout(i).time));
    dd = inn(lst,:);
    ff=f{randi(length(sampledata),1)};
    for j=1:length(ff)
        dd(:,j)=filter(1,ff{j},dd(:,j));
    end
    dd=dd+ones(size(dd,1),1)*mu{randi(length(sampledata),1)};

    X=dataout(i).getStimMatrix;
    beta=mvnrnd(sampleStats(randi(length(sampledata),1)).beta,sampleStats(randi(length(sampledata),1)).covb);
    dd=dd+X*beta;
    dataout(i,1).data=dd;
end


