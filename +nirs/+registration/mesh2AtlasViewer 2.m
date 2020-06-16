function mesh2AtlasViewer(mesh,outfolder)

img = mesh.convert2image;

dims=size(img);

if(~exist(fullfile(outfolder),'dir'))
    mkdir(outfolder);
end

if(~exist(fullfile(outfolder,'anatomical'),'dir'))
    mkdir(fullfile(outfolder,'anatomical'));
end
        
    
fid=fopen(fullfile(outfolder,'anatomical','headvol.vox'),'w');
fwrite(fid,img.vol(:),'uint8');
fclose(fid);

dlmwrite(fullfile(outfolder,'anatomical','headvol_dims.txt'),img.size);


writetable(table({'scalp','skull','csf','gm','wm'}),...
    fullfile(outfolder,'anatomical','headvol_tiss_type.txt'),...
    'FileType','text','WriteVariableName',false);


RAS=[eye(3) -img.origin'; 0 0 0 1];
dlmwrite(fullfile(outfolder,'anatomical','headvol2ras.txt'),RAS);

tbl=mesh.fiducials;

writetable(table(tbl.Name),...
    fullfile(outfolder,'anatomical','refpts_labels.txt'),...
    'FileType','text','WriteVariableName',false)

writetable(table(tbl.X,tbl.Y,tbl.Z),...
    fullfile(outfolder,'anatomical','refpts.txt'),...
    'FileType','text','WriteVariableName',false)



Tform=[diag(img.dim) img.origin'; 0 0 0 1];
dlmwrite(fullfile(outfolder,'anatomical','pialsurf2vol.txt'),Tform);
dlmwrite(fullfile(outfolder,'anatomical','refpts2vol.txt'),Tform);
dlmwrite(fullfile(outfolder,'anatomical','headsurf2vol.txt'),Tform);


write_surf(fullfile(outfolder,'anatomical','headsurf.mesh'),mesh(1).nodes,mesh(1).faces+1);  % add 1 because FS will strip off 1 and AtlasViewer needs to start at one (hack around freesurfer)
write_surf(fullfile(outfolder,'anatomical','pialsurf.mesh'),mesh(end).nodes,mesh(end).faces+1);  % add 1 because FS will strip off 1 and AtlasViewer needs to start at one (hack around freesurfer)

% add the Atlas labels (if avaliable);
aal_ll={};
aal_fv={};
cnt=1;
for i=1:mesh(end).labels.count
    if(~strcmp(mesh(end).labels.keys{i},'Structures'))
        
        a=mesh(end).labels(mesh(end).labels.keys{i});
        for j=1:length(a.Label)
            
            f=mesh(end).faces(find(all(ismember(mesh(end).faces,a.VertexIndex{j}),2)),:);
            [~,f]=ismember(f,a.VertexIndex{j});
            aal_ll{cnt}=a.Label{j};
            aal_fv{cnt}=struct('vertices',mesh(end).nodes(a.VertexIndex{j},:),'faces',f);
            cnt=cnt+1;
        end
    end
end

save(fullfile(outfolder,'anatomical','labelssurf.mat'),'aal_fv','aal_ll');
dlmwrite(fullfile(outfolder,'anatomical','labelssurf2vol.txt'),Tform);






