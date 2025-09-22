function probe2JSON(probe,filename,filename2)
% 
% if(isa(probe,'nirs.core.Probe'))
%     % this is a 2D probe and nothing to do
%     return
% end

link=probe.link;
flds=link.Properties.VariableNames;
ll=struct; 
ord={'name','type','source','detector','wavelength_nominal','units'};
for i=1:length(ord); 
    ll=setfield(ll,ord{i},link.(ord{i})); 
end;
flds={flds{find(~ismember(flds,ord))}};
for i=1:length(flds); 
    ll=setfield(ll,flds{i},link.(flds{i})); 
end;
ll.type=repmat({'NIRSCWAMPLITUDE'},height(link),1);
link=struct2table(ll);

if(strcmp(class(probe),'nirs.core.Probe1020'))
    
    
    writetable(link,[filename2 '_channels.tsv'],'FileType','text','Delimiter','\t');
    
    optodes=probe.optodes;
    optodes.name=optodes.Name; optodes.Name=[];
    optodes.type=optodes.Type; optodes.Type=[];
    optodes.x=optodes.X; optodes.X=[];
    optodes.y=optodes.Y; optodes.Y=[];
    optodes.z=optodes.Z; optodes.Z=[];
    
    optodesR=probe.optodes_registered;
    optodes.template_x=optodesR.X; 
    optodes.template_y=optodesR.Y; 
    optodes.template_z=optodesR.Z; 
    
    optodes.units=optodes.Units; optodes.Units=[];
   
    
    writetable(optodes,[filename '_optodes.tsv'],'FileType','text','Delimiter','\t');
    flds=optodes.Properties.VariableNames;
    fid=fopen([filename '_optodes.json'],'w');
    optodes.Properties.VariableDescriptions=repmat({'unknown'},length(flds),1);
    optodes.Properties.VariableDescriptions{find(ismember(flds,'type'))}='source or detector';
    optodes.Properties.VariableDescriptions{find(ismember(flds,'units'))}='spatial units';
    optodes.Properties.VariableDescriptions{find(ismember(flds,'x'))}='X position in 2D';
    optodes.Properties.VariableDescriptions{find(ismember(flds,'y'))}='Y position in 2D';
    optodes.Properties.VariableDescriptions{find(ismember(flds,'z'))}='Z position in 2D';
    optodes.Properties.VariableDescriptions{find(ismember(flds,'template_x'))}='X position in 3D';
    optodes.Properties.VariableDescriptions{find(ismember(flds,'template_y'))}='Y position in 3D';
    optodes.Properties.VariableDescriptions{find(ismember(flds,'template_z'))}='Z position in 3D';

    fprintf(fid,'{\n');
    for ii=1:length(flds)
        fprintf(fid,'\t"%s":{\n',flds{ii});
        fprintf(fid,'\t\t"description": "%s",\n',optodes.Properties.VariableDescriptions{ii});
    end
    fprintf(fid,'}');
    fclose(fid);

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
    fprintf(fid,'\t"IntendedFor": "%s"\n',['anat/T1w.nii.gz']);
    fprintf(fid,'}');
    fclose(fid);
elseif(strcmp(class(probe),'eeg.core.Probe'))
    
    electrodes=probe.electrodes;
    electrodes.X=[];
    electrodes.Y=[];
    electrodes.Z=[];
    electrodes.Units=[];
    electrodes.name=electrodes.Name; electrodes.Name=[]; 
    electrodes.type=repmat({'EEG'},height(electrodes),1);
    electrodes.units=repmat({'scaled'},height(electrodes),1);
    electrodes.status=repmat({'unknown'},height(electrodes),1);
    electrodes.status_description=repmat({' '},height(electrodes),1);
    
    writetable(electrodes,[filename '_channels.tsv'],'FileType','text','Delimiter','\t');   
    
    
    electrodes=probe.electrodes;
    electrodes.x=electrodes.X; electrodes.X=[];
    electrodes.y=electrodes.Y; electrodes.Y=[];
    electrodes.z=electrodes.Z; electrodes.Z=[];
    electrodes.units=electrodes.Units; electrodes.Units=[];
    writetable(electrodes,[filename '_electrodes.tsv'],'FileType','text','Delimiter','\t');
    
   
else
    writetable(link,[filename2 '_channels.tsv'],'FileType','text','Delimiter','\t');
    
    % 
    % optodes=probe.optodes;
    % optodes.name=optodes.Name; optodes.Name=[];
    % optodes.type=lower(optodes.Type); optodes.Type=[];
    % optodes.x=optodes.X; optodes.X=[];
    % optodes.y=optodes.Y; optodes.Y=[];
    % optodes.z=optodes.Z; optodes.Z=[];
    % optodes.units=optodes.Units; optodes.Units=[];
    % writetable(optodes,[filename '_optodes.tsv'],'FileType','text','Delimiter','\t');
    % 

 
end


return