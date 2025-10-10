function id=get_subjectID(data)
% This searches the demographics for possible ways to define the subject ID
% and returns 


if(length(data)>1)
    for i=1:length(data)
        id{i,1}=nirs.reports.helper.get_subjectID(data(i));
    end
    return
end

names={'id','subjid','subject','name','uuid'};

[a,b]=ismember(lower(data(1).demographics.keys),names);

if(any(a))
    id=data(1).demographics.values{min(find(b>0))};
else
    id='?';
end
