function nirsite2atlasviewer(montage_folder)
% montage_folder: path of montage folder created with NIRx NIRSite 2.0
%
% nirsite2atlasviewer reads files containing the channel list ('*Channels.txt') and optodes coordinates ('*Optodes.txt') in the NIRx brain space
% and outputs files 'digpts.txt' and 'probe.SD' required by AtlasViewer for computing
% forward model and sensitivity map. The output files are saved in the same
% montage folder.

%% Create AtlasViewer-compatible digitized points (digpts.txt)

% Loading NIRSite 2.0 montage coordinates

file = dir(fullfile([montage_folder filesep '*Optodes.txt']));
optodes = importdata([file.folder filesep file.name]);
coords = optodes.data;
labels = optodes.textdata;

ns = 0;
nd = 0;
for i = 1:length(labels)
    if ismember('S',char(labels(i)))
        ns = ns+1;
    end
    if ismember('D',char(labels(i)))
        nd = nd+1;
    end
    
end

% Create digitized points file (digpts.txt)

fid=fopen([montage_folder filesep 'digpts.txt'],'w');
fprintf(fid,'nz: 0.4 85.9 -47.6\n');
fprintf(fid,'ar: 83.9 -16.6 -56.7\n');
fprintf(fid,'al: -83.8 -18.6 -57.2\n');
fprintf(fid,'cz: -0.461 -8.416 101.365\n');
fprintf(fid,'iz: 0.2 -120.5 -25.8\n');
for i = 1:ns
    fprintf(fid,'s%d: %.1f %.1f %.1f\n',i,coords(i,1),coords(i,2),coords(i,3));
end
for i = 1:nd
    fprintf(fid,'d%d: %.1f %.1f %.1f\n',i,coords(ns+i,1),coords(ns+i,2),coords(ns+i,3));
end
fclose(fid);


%% Create AtlasViewer-compatible probe condiguration (probe.SD)

% Loading NIRSite montage coordinates

file = dir(fullfile([montage_folder filesep '*Channels.txt']));
channels = importdata([file.folder filesep file.name]);

MeasList = [channels(:,2:3) ones(length(channels),2)]; 
MeasList2 = MeasList;
MeasList2(:,4) = 2;
SD.MeasList = [MeasList; MeasList2];

% convert for NIRS toolbox format
SD.SrcPos = coords(1:ns,:);
SD.DetPos = coords(ns+1:ns+nd,:);
SD.nSrcs = ns;
SD.nDets = nd;
save([montage_folder filesep 'probe.SD'],'SD')
