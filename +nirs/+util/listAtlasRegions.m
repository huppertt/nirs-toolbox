function lst=listAtlasRegions
%This function returns a list of the regions in the AAL atlas

aalLabels=load(which('ROI_MNI_V5_List.mat'));
lst=strvcat(aalLabels.ROI.Nom_L);

return
