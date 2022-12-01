function data = remove_demographic(data,fld)

if(ischar(fld))
    fld={fld};
end

demo=nirs.createDemographicsTable(data);
demo(:,find(ismember(demo.Properties.VariableNames,fld)))=[];

fld=demo.Properties.VariableNames;

for i=1:length(data)
    data(i).demographics=Dictionary;
    for j=1:length(fld)
        data(i).demographics(fld{j})=demo.(fld{j})(i);
    end
end