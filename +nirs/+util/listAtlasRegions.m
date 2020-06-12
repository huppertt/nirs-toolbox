function lst=listAtlasRegions(labels)
%This function returns a list of the regions in the AAL atlas

if(nargin==0)
    aalLabels=load(which('ROI_MNI_V5_List.mat'));
    lst=strvcat(aalLabels.ROI.Nom_L);
else
    keys=labels.keys;
    lst={};
    for i=1:length(keys)
        l=labels(keys{i});
        lst={lst{:} l.Label{:}};
    end
    lst=strvcat(lst);
end
return
