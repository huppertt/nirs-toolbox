function data=clean_up_demographics(data)

tbl=nirs.createDemographicsTable(data);

flds=tbl.Properties.VariableNames;

for i=1:length(flds)
    if(iscell(tbl.(flds{i})))
        for j=1:height(tbl)
            if(strcmp(tbl.(flds{i}){j},'null'))
                tbl.(flds{i}){j}=NaN;
            end
        end
        try;
            if(isnumeric(cell2mat(tbl.(flds{i}))))
                tbl.(flds{i})=cell2mat(tbl.(flds{i}));
            end
        end
    elseif(isstruct(tbl.(flds{i})))
        tmp={};
        for j=1:height(tbl)
            tmp{j,1}=jsonencode(tbl.(flds{i})(j),'PrettyPrint',false);
        end
        tbl.(flds{i})=tmp;

    end
end

for i=1:length(data)
    for j=1:length(flds)
        data(i).demographics(flds{j})=tbl(i,:).(flds{j});
    end
end
