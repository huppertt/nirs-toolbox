function DataSet2BIDS(data,folder,name,LocalBIDSinfo,ReadMe,deID)

if(nargin<4 || isempty(LocalBIDSinfo))
    LocalBIDSinfo=struct;
end

if(nargin<5 || isempty(ReadMe))
    ReadME='Raw data from NIRS toolbox';
end

LocalBIDSinfo.BIDSVersion='1.8.0';
LocalBIDSinfo.License='CC0';
LocalBIDSinfo.DatasetType='raw';
LocalBIDSinfo.GeneratedBy.Name='nirs-toolbox';
LocalBIDSinfo.GeneratedBy.Version='2.0.0';

 



if(~exist(folder,'dir'))
    mkdir(folder);
end

if(nargin<3)
    name=folder;
end

if(nargin<5)
    deID=false;
end

data = nirs.util.addUUID(data);

if(deID)
    data=nirs.util.deidentify(data);
end



%% Make the dataset_description.json file
fid=fopen(fullfile(folder,'dataset_description.json'),'w');
fprintf(fid,'{\n');
fprintf(fid,'  "Name": "%s",\n',name);
fprintf(fid,'  "DatasetDOI": "%s",\n',datestr(now));
flds=fields(LocalBIDSinfo);
for f=1:length(flds);
    if(iscellstr(LocalBIDSinfo.(flds{f})))
            fprintf(fid,'  "%s": [\n',flds{f});
            for j=1:length(LocalBIDSinfo.(flds{f}))
                 fprintf(fid,'  \t "%s"',LocalBIDSinfo.(flds{f}){j});
                if(j~=length(LocalBIDSinfo.(flds{f})))
                    fprintf(fid,',\n');
                end
            end
            if(f~=length(flds))
                fprintf(fid,'\n  ],\n');
            else
                fprintf(fid,'\n  ]\n');
            end
    elseif(isstruct(LocalBIDSinfo.(flds{f})))
        fprintf(fid,' "%s": [\n',flds{f});
        stct=LocalBIDSinfo.(flds{f});
        flds2=fields(stct);
        fprintf(fid,' \t\t{\n');
        for f2=1:length(flds2)
            if(f2~=length(flds2))
                fprintf(fid,'  \t\t\t"%s": "%s",\n',flds2{f2},stct.(flds2{f2}));
            else
                fprintf(fid,'  \t\t\t"%s": "%s"\n',flds2{f2},stct.(flds2{f2}));
            end
        end
        fprintf(fid,' \t\t}\n');
         if(f~=length(flds))
            fprintf(fid,'  \t],\n');
         else
            fprintf(fid,'  \t]\n');
         end
    else
        if(f~=length(flds))
            fprintf(fid,'  "%s": "%s",\n',flds{f},LocalBIDSinfo.(flds{f}));
        else
            fprintf(fid,'  "%s": "%s"\n',flds{f},LocalBIDSinfo.(flds{f}));
        end
    end
end
fprintf(fid,'}');
fclose(fid);

fid=fopen(fullfile(folder,'README.md'),'w');
fprintf(fid,'%s',ReadME);
fclose(fid);


%% Make the participants files
demo = nirs.createDemographicsTable(data);

data = nirs.bids.add_missing_BIDS_info(data);


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
    Task=lower(Task);
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
    Session=lower(Session);
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

[~,lstI,lstL]=unique(demo.participant_id);
demo2=demo(lstI,stable);


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
    folder2=fullfile(folder,demo2.participant_id{i});
    
    if(~exist(folder2,'dir'))
       system(['mkdir -p ' folder2]);
    end
    lst2=find(lstL==i);   
    
    filenamesAll={};
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
        if(iscellstr(run))
            run=run{1};
        end
        run=['0' run];
        run=run(end-1:end);
        run={run};
        
        folder3=folder2;
        if(length(unique({Session{lst2}}))>1)
            folder3=fullfile(folder3,['ses-' Session{lst2(j)}]);
            lst3=find(ismember({Session{lst2}},Session{lst2(j)}));
            this_sess=Session{lst2(j)};
        else
            this_sess=[];
        end
        
        system(['mkdir -p ' folder3]);

        filenamestmp=Data2BIDS(data(lst2(j)),folder3,demo2.participant_id{i},lower(this_sess),lower(Task{lst2(j)}),run);
        filenamesAll{j,1}=filenamestmp{1};
    end
    demo3=demo(lst2,:);
    filenamesAll=strrep(filenamesAll,[folder filesep],'');
    filenamesAll=strrep(filenamesAll,[demo2.participant_id{i} filesep],'');
    
    demo3=[table(filenamesAll,'VariableNames',{'filename'}) demo3];
    stable=true(length(demo3.Properties.VariableNames),1);
    for ii=1:length(stable)

        try
            if(~all(isnan(table2array(demo3(:,ii)))))
                stable(ii)=stable(ii) & (height(unique(demo3(:,ii)))==1);
            end
        catch
            stable(ii)=stable(ii) & (height(unique(demo3(:,ii)))==1);
        end
    end
    if(~all(stable))
        demo3=demo3(:,~stable);

        fid=fopen(fullfile(folder2,[demo2.participant_id{i} '_scans.json']),'w');
        fprintf(fid,'{\n');
        for ii=1:length(demo3.Properties.VariableNames)
            if(ii>1)
                fprintf(fid,',\n');
            end
            if(isempty(demo3.Properties.VariableDescriptions))
                demo3.Properties.VariableDescriptions=repmat({' '},length(demo3.Properties.VariableNames),1);
            end
            if(isempty(demo3.Properties.VariableUnits))
                demo3.Properties.VariableUnits=repmat({' '},length(demo3.Properties.VariableNames),1);
            end
            fprintf(fid,'\t"%s": {\n\t\t"Description": "%s",\n\t\t"Units": "%s"}\n\t',...
                demo3.Properties.VariableNames{ii},...
                demo3.Properties.VariableDescriptions{ii},...
                demo3.Properties.VariableUnits{ii});
        end
        fprintf(fid,'\n}');
        fclose(fid);

        writetable(demo3,fullfile(folder2,[demo2.participant_id{i} '_scans.tsv']),'FileType','text','Delimiter','\t');
    end
end




    






