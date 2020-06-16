function probe2JSON(probe,filename)
% 
% if(isa(probe,'nirs.core.Probe'))
%     % this is a 2D probe and nothing to do
%     return
% end

writetable(probe.optodes,[filename '_channels.tsv'],'FileType','text','Delimiter','\t');

if(strcmp(class(probe),'nirs.core.Probe1020'))
    writetable(probe.optodes_registered,[filename '_optodes.tsv'],'FileType','text','Delimiter','\t');
    
    mesh=probe.getmesh;
    fidc=mesh(1).fiducials;
    
    fid=fopen([filename '_coordsystem.json'],'w');
    
    fprintf(fid,'{\n');
    fprintf(fid,'\t"NIRSCoordinateSystem": "%s",\n','MNI');
    fprintf(fid,'\t"NIRSCoordinateUnits": "%s",\n','mm');
    fprintf(fid,'\t"AnatomicalLandmarkCoordinates": {\n');
    for i=1:height(fidc)-1
        fprintf(fid,'\t\t"%s": [%f,\t%f,\t%f],\n',fidc.Name{i},fidc.X(i),fidc.Y(i),fidc.Z(i));
    end
    fprintf(fid,'\t\t"%s": [%f,\t%f,\t%f]\n',fidc.Name{end},fidc.X(end),fidc.Y(end),fidc.Z(end));
    fprintf(fid,'\t},\n');
    fprintf(fid,'\t"AnatomicalLandmarkCoordinateSystem": "%s",\n','T1w');
    fprintf(fid,'\t"AnatomicalLandmarkCoordinateUnits": "%s",\n','mm');
    fprintf(fid,'\t"IntendedFor": "%s",\n',['anat/T1w.nii.gz']);
    fprintf(fid,'}');
    fclose(fid);
    
    
end


return