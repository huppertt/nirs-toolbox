function DataSet2BIDS(data,folder,name,deID)
BIDSver='pre-release 0.1';

if(~exist(folder,'dir'))
    mkdir(folder);
end

if(nargin<3)
    name=folder;
end

if(nargin<4)
    deID=false;
end

data = nirs.util.addUUID(data);

if(deID)
    data=nirs.util.deidentify(data);
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


if(isempty(demo))
    for i=1:length(data)
        if(isempty(data(i).description))
            data(i).description=' ';
        end
        data(i).demographics('description')=data(i).description;
    end
    demo = nirs.createDemographicsTable(data);
end

lst=find(ismember(lower(demo.Properties.VariableNames),{'task','tasks','scans'}));
if(~isempty(lst))
    for i=1:height(demo)
        Task{i,1}=[];
        for j=1:length(lst)
            try
                Task{i}=[Task{i} demo.(demo.Properties.VariableNames{lst(j)}){i}];
            end
        end
        
    end
    demo(:,lst)=[];
else
    for i=1:height(demo)
        Task{i,1}=[];
    end
end
lst=find(ismember(lower(demo.Properties.VariableNames),{'session','sessions','study'}));
if(~isempty(lst))
    for i=1:height(demo)
        Session{i,1}=[];
        for j=1:length(lst)
            try
                Session{i}=[Session{i} demo.(demo.Properties.VariableNames{lst(j)}){i}];
            end
        end
        
    end
    demo(:,lst)=[];
else
    for i=1:height(demo)
        Session{i,1}=[];
    end
end
% 
 demo.UUID=[];
% demo.SubjID=[];
% demo.MeasurementTime=[];
% demo.MeasurementDate=[];


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
    fprintf(fid,'\t"%s": {\n\t\t"Description": "%s",\n\t\t"Units": "%s"}\n\t',...
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

    if(isa(data,'nirs.core.Data'))
        folder2=fullfile(folder,demo2.participant_id{i},'nirs');
    elseif(isa(data,'eeg.core.Data'))
        folder2=fullfile(folder,demo2.participant_id{i},'eeg');
    else
        error('unknown BIDS type');
    end
    
    if(~exist(folder2,'dir'))
        mkdir(fullfile(folder,demo2.participant_id{i}));
        mkdir(folder2);
    end
    lst2=find(lst==i);
    for j=1:length(lst2)
        
        
        if(length(lst2)>1)
            run=num2str(j);
        else
            run=[];
        end
               
        
        Data2BIDS(data(lst2(j)),folder2,demo2.participant_id{i},Session{lst2(j)},Task{lst2(j)},run);
    end
end




    






