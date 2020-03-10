function data2JSON(data,filename,task)

fid=fopen(filename,'w');
fprintf(fid,'{\n');

snirf=nirs.util.validateSNIRF(data);

fprintf(fid,'\t"TaskName": %s\n',task);
fprintf(fid,'\t"Version": %s\n',snirf.formatVersion);
fprintf(fid,'\t"MeasurementDate": %s\n',snirf.nirs.metaDataTags.MeasurementDate);
fprintf(fid,'\t"MeasurementTime": %s\n',snirf.nirs.metaDataTags.MeasurementTime);
fprintf(fid,'\t"SNIRF_createDate": %s\n',snirf.nirs.metaDataTags.SNIRF_createDate);
fprintf(fid,'\t"SNIRF_createTime": %s\n',snirf.nirs.metaDataTags.SNIRF_createTime);
fprintf(fid,'\t"filedescription": %s\n',snirf.nirs.metaDataTags.filedescription);
fprintf(fid,'}');
fclose(fid);

return