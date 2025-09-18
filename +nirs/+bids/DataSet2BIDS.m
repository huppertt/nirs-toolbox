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



lst=find(ismember(lower(demo.Properties.VariableNames),{'task','tasks','cond'}));
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
lst=find(ismember(lower(demo.Properties.VariableNames),{'session','sessions','study','visit','timepoint'}));
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

lst=find(ismember(lower(demo.Properties.VariableNames),{'subject','subjid','id','subjectid'}));

if(~isempty(lst))
    for i=1:height(demo)
        participant_id{i,1}=['sub-' demo.(demo.Properties.VariableNames{lst(1)}){i}];
    end
else
    for i=1:height(demo)
        str=['000' num2str(i)];
        participant_id{i,1}=['sub-' str(end-4:end)];
    end
end
demo=[table(participant_id) demo];

upI=unique(demo.participant_id);
stable=true(length(demo.Properties.VariableNames),1);
for i=1:length(stable)
    for j=1:length(upI);
        lst=find(ismember(demo.participant_id,upI{j}));
        try
        if(~all(isnan(table2array(demo(lst,i)))))
            stable(i)=stable(i) & (height(unique(demo(lst,i)))==1);
        end
        catch
            stable(i)=stable(i) & (height(unique(demo(lst,i)))==1);
        end
    end
end

[demo2,~,lstL]=unique(demo(:,stable));


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
    lst2=find(lstL==i);   

    for j=1:length(lst2)
        
        lst=find(ismember(lower(demo.Properties.VariableNames),{'scans','scan','run'}));
        if(~isempty(lst))
            run=strtrim(demo.(demo.Properties.VariableNames{lst(1)})(lst2(j)));
        else
            if(length(lst2)>1)
                run=num2str(j);
            else
                run=[];
            end
        end   

        folder3=folder2;
        if(length(unique({Session{lst2}}))>1)
            folder3=fullfile(folder3,['ses-' Session{lst2(j)}]);
            lst3=find(ismember({Session{lst2}},Session{lst2(j)}));
        else
            lst3=1:length(lst3);
        end
        if(length(unique({Task{lst2(lst3)}}))>1)
            folder3=fullfile(folder3,['task-' Task{lst2(j)}]);
        end
        system(['mkdir -p ' folder3]);

        Data2BIDS(data(lst2(j)),folder3,demo2.participant_id{i},Session{lst2(j)},Task{lst2(j)},run);
    end
    demo3=demo(lst2,:);
    stable=true(length(demo3.Properties.VariableNames),1);
    for i=1:length(stable)

        try
            if(~all(isnan(table2array(demo3(:,i)))))
                stable(i)=stable(i) & (height(unique(demo3(:,i)))==1);
            end
        catch
            stable(i)=stable(i) & (height(unique(demo3(:,i)))==1);
        end
    end
    if(~all(stable))
        demo3=demo3(:,~stable);

        fid=fopen(fullfile(folder2,'scans.json'),'w');
        fprintf(fid,'{\n');
        for i=1:length(demo3.Properties.VariableNames)
            if(i>1)
                fprintf(fid,',\n');
            end
            if(isempty(demo3.Properties.VariableDescriptions))
                demo3.Properties.VariableDescriptions=repmat({' '},length(demo3.Properties.VariableNames),1);
            end
            if(isempty(demo3.Properties.VariableUnits))
                demo3.Properties.VariableUnits=repmat({' '},length(demo3.Properties.VariableNames),1);
            end
            fprintf(fid,'\t"%s": {\n\t\t"Description": "%s",\n\t\t"Units": "%s"}\n\t',...
                demo3.Properties.VariableNames{i},...
                demo3.Properties.VariableDescriptions{i},...
                demo3.Properties.VariableUnits{i});
        end
        fprintf(fid,'\n}');
        fclose(fid);

        writetable(demo3,fullfile(folder2,'scans.tsv'),'FileType','text','Delimiter','\t');
    end
end




    






