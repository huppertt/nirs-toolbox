function [] = gin_dlabels()
% ______________________________________________________________________
%
% Local Maxima Labeling : gin_dlabels.m
% Functional Maps are thresholded and the clusters and local maxima are
% extracted.
% The function is the assignment of the label of the 3D volume of interest 
% it belongs to. The three nearest anatomical regions are listed.
%
% gin_dlabels.m		B Landeau - 10/09/03 - aal for SPM2
% gin_dlabels.m		G Flandin - 21/08/15 - aal for SPM12
% ______________________________________________________________________


spm('defaults','FMRI');
try
    xSPM = evalin('base','xSPM');
    hReg = evalin('base','hReg');
catch
    [hReg,xSPM,SPM] = spm_results_ui;
    TabDat = spm_list('List',xSPM,hReg);
    assignin('base','SPM',SPM);
    assignin('base','xSPM',xSPM);
    assignin('base','hReg',hReg);
    assignin('base','TabDat',TabDat);
end
% Compute and Display Labels (and distances) for local max and for cluster
gin_list_dlabels('List',xSPM,hReg);

