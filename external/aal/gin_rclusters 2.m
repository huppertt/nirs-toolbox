function [] = gin_rclusters()
% ______________________________________________________________________
%
% Extended Local Maxima Labeling : gin_rclusters.m
% Functional Maps are thresholded and the clusters and local maxima are
% extracted.
% A spherical volume is defined around the local maxima. The radius of this
% region is chosen by the user (default value = 10 mm). For each volume of
% interest, the number of voxel overlapping this spherical region is computed.
% Note : A label "OUTSIDE" is included if some parts of the spherical regions
% fall outside any volume of interest.
%
% gin_rclusters.m		B Landeau - 10/09/03 - aal for SPM2
% ______________________________________________________________________

% get the SPM.mat file.
%[SPM,VOL,xX,xCon,xSDM] = spm_getSPM;
%[SPM,xSPM] = spm_getSPM;
try
    xSPM = evalin('base','xSPM');
catch
    [SPM,xSPM] = spm_getSPM;
end

% Compute Labels (in %) for local max and for cluster
% from a spherical region defined by the user.
% and Display the Results Table :
% 		local maxima		label		%
gin_list_plabels('List',xSPM);
