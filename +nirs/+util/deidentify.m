function data = deidentify(data,keep)

if(nargin<2)
    keep={'uuid','age','gender','task','session' };
end


data=nirs.util.addUUID(data);

demographics=nirs.createDemographicsTable(data);

genderAlias={'gender','sex'};
lst=find(ismember(lower(demographics.Properties.VariableNames),genderAlias));
if(~isempty(lst))
for i=1:height(demographics)
    g='unknown';
    for j=1:length(lst)
        if(ismember(lower(demographics.(demographics.Properties.VariableNames{lst(j)})(i)),{'male','m'}))
            g='male';
        elseif(ismember(lower(demographics.(demographics.Properties.VariableNames{lst(j)})(i)),{'female','f'}))
            g='female';
        elseif(ismember(lower(demographics.(demographics.Properties.VariableNames{lst(j)})(i)),{'nonbinay','non-binay','none','cis','trans'}))
            g='other';
        end
    end
    gender{i,1}=g;
    
end
demographics(:,lst)=[];
demographics=[demographics table(gender)];
end

lst=find(ismember(lower(demographics.Properties.VariableNames),'age'));
if(~isempty(lst))
    for i=1:height(demographics)
        a=demographics.(demographics.Properties.VariableNames{lst})(i);
        if(iscellstr(a))
            a=str2num(a{1});
        end
        if(isempty(a))
            a=NaN;
        end
        age(i,1)=a;
    end
    demographics(:,lst)=[];
    demographics=[demographics table(age)];
end

% remove anything in the removeFlds list
lst=find(~ismember(lower(demographics.Properties.VariableNames),lower(keep)));
demographics(:,lst)=[];

flds=demographics.Properties.VariableNames;
for i=1:height(demographics)
    data(i).demographics=Dictionary;
    for j=1:length(flds);
        
        data(i).demographics(flds{j})=demographics.(flds{j})(i);
    end
    data(i).description='deidentified';
end



return