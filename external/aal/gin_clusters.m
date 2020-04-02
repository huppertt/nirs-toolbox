function [] = gin_clusters()
% ______________________________________________________________________
%
% Cluster Labeling : gin_clusters.m
% Functional Maps are thresholded and the clusters and local maxima are
% extracted.
% The same procedure as in extented local maxima labelleing is used, 
% replacing the spherical region by the functional cluster.
%
% gin_clusters.m		B Landeau - 09/09/03 - aal for SPM2
% ______________________________________________________________________

% get the SPM.mat file.
try
    xSPM = evalin('base','xSPM');
catch
    [SPM,xSPM] = spm_getSPM;
end

% Compute Labels (in %) for cluster
% from a spherical region defined by the user.
% and Display the Results Table :
% 		local maxima		label		%
gin_clusters_plabels('List',xSPM);
