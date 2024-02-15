function data2JSON_EEG(data,filename,task)

fid=fopen(filename,'w');
fprintf(fid,'{\n');

file=data.description;
file=[strtok(file,'.') '.vhdr'];
[hdr.fs, hdr.label,hdr.meta]=bva_readheader(file);


fprintf(fid,'\t"TaskName": "%s",\n',task);
fprintf(fid,'\t"SamplingFrequency": %f,\n',hdr.fs);
fprintf(fid,'\t"SoftwareFilters": "%s",\n',"n/a");
fprintf(fid,'\t"EEGChannelCount": %d,\n',height(data.probe.link));
fprintf(fid,'\t"EOGChannelCount": %d,\n',length(find(ismember(lower(hdr.label),'eog'))));
fprintf(fid,'\t"ECGChannelCount": %d\n',length(find(ismember(lower(hdr.label),'ecg'))));
%fprintf(fid,'\t"EEGReference": %s\n',snirf.nirs.metaDataTags.SNIRF_createTime);
%fprintf(fid,'\t"PowerLineFrequency": %s\n',snirf.nirs.metaDataTags.SNIRF_createTime);
%fprintf(fid,'\t"filedescription": %s\n',snirf.nirs.metaDataTags.filedescription);
fprintf(fid,'}');
fclose(fid);

return