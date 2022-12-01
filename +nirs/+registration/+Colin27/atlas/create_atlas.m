function atlas =create_atlas(file,version1)

nii=load_nii(file);


T1=[ nii.hdr.hist.srow_x;  nii.hdr.hist.srow_y;  nii.hdr.hist.srow_z; 0 0 0 1];
T =[ 0.9964    0.0178    0.0173   -0.0000
-0.0169    0.9957   -0.0444   -0.0000
-0.0151    0.0429    1.0215    0.0000
-0.4232  -17.5022   11.6967    1.0000];  %MNI to talriach

if(nargin==2 && version1)
    mesh=nirs.registration.Colin27.mesh;
else
    % version 2 is in the wavelet space
    mesh=nirs.registration.Colin27.mesh_V2;
end

n=mesh(end).nodes;
n(:,4)=1;
n1=n*inv(T)*inv(T1');
n1(:,4)=[];

lst=find(nii.img~=0);
[n2(:,1),n2(:,2),n2(:,3)]=ind2sub(size(nii.img),lst);

[k,d]=dsearchn(n2,n1);


lst2=unique(nii.img(lst(k)));
for i=1:length(lst2); 
    atlas.VertexIndex{i,1}=find(nii.img(lst(k))==lst2(i)); 
end;

[p,file]=fileparts(file);
[~,file]=fileparts(file);

tbl=readtable([p file '.txt']);
tbl=tbl(lst2,:);
atlas.Label=tbl.Var2;
atlas.Region=tbl.Var2;



aal=load('ROI_MNI_V5_Border_modified.mat');
%aalLabels=load(which('ROI_MNI_V5_List.mat'));
aal.BORDER_XYZ(1,:)=aal.BORDER_XYZ(1,:)*2-90;
aal.BORDER_XYZ(2,:)=aal.BORDER_XYZ(2,:)*2-126;
aal.BORDER_XYZ(3,:)=aal.BORDER_XYZ(3,:)*2-72;
aal.BORDER_XYZ=icbm_spm2tal(aal.BORDER_XYZ')';

aal.BORDER_XYZ(1,:)=aal.BORDER_XYZ(1,:)-2.5;
aal.BORDER_XYZ(2,:)=aal.BORDER_XYZ(2,:)+17.5;
aal.BORDER_XYZ(3,:)=aal.BORDER_XYZ(3,:)-20;

n2=aal.BORDER_XYZ';
n=mesh(end).nodes;
[k,d]=dsearchn(n2,n);


for ii=1:3
    lst2=unique(aal.BORDER_V(ii,k));
    for i=1:length(lst2)
        for j=1:length(aal.ROI{ii})
            if(aal.ROI{ii}(j).ID==lst2(i))
                atlas{ii}.VertexIndex{i,1}=find(aal.BORDER_V(ii,k)==aal.ROI{ii}(j).ID)';
                if(isempty(aal.ROI{ii}(j).Nom_L))
                    atlas{ii}.Label{i,1}=[aal.labels{ii} '_' num2str(j)]
                    atlas{ii}.Region{i,1}=[aal.labels{ii} '_' num2str(j)]
                    
                else
                    atlas{ii}.Label{i,1}=aal.ROI{ii}(j).Nom_L;
                    atlas{ii}.Region{i,1}=aal.ROI{ii}(j).Nom_C;
                end
            end
        end
    end;
end
fwdBEM.mesh(end).labels('aal')=atlas{1}
fwdBEM.mesh(end).labels('Brodmann (MRIcron)')=atlas{2};
fwdBEM.mesh(end).labels('Brodmann (Talairach daemon)')=atlas{3};
