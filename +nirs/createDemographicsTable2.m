function tbl = createDemographicsTable2( data )
% This is a simplified version of createDemographicsTable that is around 40x faster

% Get list of all fields in first file
if iscell(data(1).demographics)
    fields = data(1).demographics{1}.keys(:)';
    for j = 2:length(data(1).demographics)
        fields = [fields strcat(data(1).demographics{1}.keys(:)','_',num2str(j))];
    end
else
    fields = data(1).demographics.keys(:)';
end
[fields,sidx] = sort(fields);
[~,sidx] = sort(sidx);

% Extract data to cell matrix
outcell = cell(length(data),length(fields));
for i=1:length(data)
    
    if iscell(data(i).demographics)
        keys = data(i).demographics{1}.keys(:)';
        values = data(i).demographics{1}.values;
        for j = 2:length(data(i).demographics)
            keys = [keys strcat(data(i).demographics{j}.keys(:)','_',num2str(j))];
            values = [values data(i).demographics{j}.values];
        end
    else
        keys = data(i).demographics.keys(:)';
        values = data(i).demographics.values;
    end
    
    [keys,tmpidx] = sort(keys);
    assert(isequal(fields,keys),'Inconsistent keys');
        
    outcell(i,:) = values(tmpidx);
    
end

% Convert to table
tbl=cell2table(outcell(:,sidx),'VariableNames',fields(sidx));

% Prune duplicate hyperscanning fields
if iscell(data(1).demographics) && length(data(1).demographics)>1
    origfields = sort(data(1).demographics{1}.keys);
    for i = 1:length(origfields)
        
        % Check if fields are duplicates
        isduplicate = true;
        for j = 2:length(data(1).demographics)
            tmpkey = strcat(origfields{i},'_',num2str(j));
            if ~isequaln(tbl.(origfields{i}),tbl.(tmpkey))
                isduplicate = false;
            end
        end
        
        % Remove duplicate fields
        if isduplicate
            for j = 2:length(data(1).demographics)
                tmpkey = strcat(origfields{i},'_',num2str(j));
                tbl.(tmpkey) = [];
            end 
        end
        
    end
end

end