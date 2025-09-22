function stim2JSON(stimdict,filename)

if(~isa(stimdict,'Dictionary'))
    tmp=Dictionary;
    tmp(stimdict.name)=stimdict;
    stimdict=tmp;
    clear tmp;
end

s=struct;
s.onset=[];
s.duration=[];
s.amplitude=[];
s.name={};

for i=1:stimdict.count
    st=stimdict(stimdict.keys{i});
    for j=1:length(st.onset)
        s.onset(end+1,1)=st.onset(j);
        s.duration(end+1,1)=st.dur(j);
        s.amplitude(end+1,1)=st.amp(j);
        s.name{end+1,1}=st.name;
        
    end
end

s=struct2table(s);
s=sortrows(s);

fid=fopen(filename,'w');
fprintf(fid,'{\n');
% fprintf(fid,'\t"onset":{\n');
% fprintf(fid,'\t\t"description": "onset of event",\n');
% fprintf(fid,'\t\t"units": "seconds"},\n');
% fprintf(fid,'\t"duration":{\n');
% fprintf(fid,'\t\t"description": "duration of event",\n');
% fprintf(fid,'\t\t"units": "seconds"},\n');
fprintf(fid,'\t"amplitude":{\n');
fprintf(fid,'\t\t"description": "amplitude of event"},\n');
fprintf(fid,'\t"name": {\n');
fprintf(fid,'\t\t"description": "name of task"},\n');
fprintf(fid,'\t"StimulusPresentation": {\n');
fprintf(fid,'\t\t"SoftwareName": "unknown",\n');
fprintf(fid,'\t\t"SoftwareVersion": "unknown"}\n');
fprintf(fid,'}');

% To prevent writetable from sometimes using numbers with terrible
% precision (e,g 0.9999999999999999).  Change to a string for the write
s.onset=strtrim(cellstr(num2str(s.onset)));
s.duration=strtrim(cellstr(num2str(s.duration)));
s.amplitude=strtrim(cellstr(num2str(s.amplitude)));

writetable(s,[filename(1:strfind(filename,'.json')-1) '.tsv'],...
    'FileType','text','Delimiter','\t');
