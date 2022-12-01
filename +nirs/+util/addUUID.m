function data = addUUID(data)

% List of possible ways to state the name
subjNamesAlias={'subject','subjid','name','id'};

demographics=nirs.createDemographicsTable(data);

if(~ismember('UUID',demographics.Properties.VariableNames))
    % if the data already has a UUID, then keep that
    
    lst=find(ismember(lower(demographics.Properties.VariableNames),subjNamesAlias));
    if(~isempty(lst))
        for i=1:height(demographics)
            UUID{i,1}=[];
            for j=1:length(lst)
                if(iscellstr(demographics.(demographics.Properties.VariableNames{lst(j)})))
                    UUID{i,1}=[UUID{i,1} demographics.(demographics.Properties.VariableNames{lst(j)}){i}];
                end
            end
        end
    else
        for i=1:height(demographics)
            UUID{i,1}=num2str(i);
        end
    end
    demographics(:,lst)=[];
    demographics=[table(UUID) demographics];
    [s,id]=unique(UUID);
    
    for i=1:length(s)
        uuid = char(java.util.UUID.randomUUID);
        lst=find(ismember(UUID,s{i}));
        for j=1:length(lst)
            demographics.UUID{lst(j)}=uuid;
        end
    end
    
    for i=1:length(data)
        data(i).demographics('UUID')=demographics.UUID{i};
    end
    
end

end