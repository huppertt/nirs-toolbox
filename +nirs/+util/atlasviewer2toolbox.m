function atlasviewer2toolbox(montage_folder)
% montage_folder: path of subject or atlas folder where AtlasViewer saved
% sensitivity map
%
% atlasviewer2toolbox reads the sensitivity map created by AtlasViewer and add AnalyzIR-specific
% variables for running image reconstruction using ImageReconMFX, saving it as atlasViewer_analyzir.mat in the same subject folder.  

%% AtlasViewer -> Template atlasViewer_fwmodel\

load([montage_folder filesep 'atlasViewer.mat']);

% Adot has to be doubled too
temp = fwmodel.Adot;
fwmodel.Adot = [temp; temp];

% Post processing "sensitivity" field
sensitivity = struct();
sensitivity.name = 'sensitivity';
sensitivity.nphotons = fwmodel.nphotons;
sensitivity.timegates = fwmodel.timegates;
sensitivity.Ch = fwmodel.Ch;
sensitivity.cmThreshold = fwmodel.cmThreshold;
sensitivity.mesh = fwmodel.mesh;
sensitivity.Adot = fwmodel.Adot;

% Post processing probe links
temp = probe.mlmp;
temp2 = probe.ml;
probe.mlmp = [temp; temp];
probe.ml = [temp2; temp2];
probe.ml(length(temp2)+1:end, 4) = 2;

% remove some of the variables you defined above but don't need
clearvars temp temp2

save([montage_folder filesep 'atlasViewer_analyzir.mat'],'-mat');