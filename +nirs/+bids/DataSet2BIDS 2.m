function DataSet2BIDS(data,folder,name)
BIDSver='???';

if(~exist(folder,'dir'))
    mkdir(folder);
end

if(nargin<3)
    name=folder;
end


%% Make the dataset_description.json file
fid=fopen(fullfile(folder,'dataset_description.json'),'w');
fprintf(fid,'{\n');
fprintf(fid,'\t"Name": "%s",\n',name);
fprintf(fid,'\t"DataserDOI": "%s",\n',datestr(now));
fprintf(fid,'\t"BIDSVersion": "%s"\n',BIDSver);
fprintf(fid,'}');
fclose(fid);


%% Make the participants files
demo = nirs.createDemographicsTable(data);

[demo2,~,lst]=unique(demo);
fid=fopen(fullfile(folder,'participants.json'),'w');
fprintf(fid,'{\n');
for i=1:length(demo2.Properties.VariableNames)
    if(i>1)
        fprintf(fid,',\n');
    end
    if(isempty(demo2.Properties.VariableDescriptions))
       demo2.Properties.VariableDescriptions=repmat({' '},length(demo2.Properties.VariableNames),1);
    end
    if(isempty(demo2.Properties.VariableUnits))
       demo2.Properties.VariableUnits=repmat({' '},length(demo2.Properties.VariableNames),1);
    end
    fprintf(fid,'\t{"%s": {\n\t\t"Description": "%s",\n\t\t"Units": "%s"}\n\t}',...
        demo2.Properties.VariableNames{i},...
        demo2.Properties.VariableDescriptions{i},...
        demo2.Properties.VariableUnits{i});
end
fprintf(fid,'\n}');
fclose(fid);

for i=1:height(demo2)
    str=['000' num2str(i)];
    participant_id{i,1} =['Sub-' str(end-3:end)];
end
demo2= [table(participant_id) demo2];

writetable(demo2,fullfile(folder,'participants.tsv'),'FileType','text','Delimiter','\t');

for i=1:height(demo2)
    folder2=fullfile(folder,demo2.participant_id{i});
    if(~exist(folder2,'dir'))
        mkdir(folder2);
    end
    nirs.bids.Data2BIDS(data(find(lst==i)),folder2,demo2.participant_id{i});
end




    






