function lst=listAtlasRegions(labels)
%This function returns a list of the regions in the AAL atlas

if(nargin==0)
    mesh=nirs.registration.Colin27.mesh_V2;
    lst=nirs.util.listAtlasRegions(mesh(end).labels);
    %aalLabels=load(which('ROI_MNI_V5_List.mat'));
    %lst=strvcat(aalLabels.ROI.Nom_L);
else
    keys=labels.keys;
    lst={};
    for i=1:length(keys)
        l=labels(keys{i});
        for j=1:length(l.Label)
            if(~isempty(l.Label{j}))
            lst{end+1}=l.Label{j};
            end
        end
    end
    lst=strvcat(lst);
end
return
