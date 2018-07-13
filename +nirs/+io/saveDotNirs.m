function saveDotNirs( data , out_dir , filename_pattern )
%% Export data to .nirs files
% nirs.io.saveDotNirs( Data , [ output_directory , filename_pattern ] )
%
% Filename pattern: File index or demographic fields delimited by %variable%
%
% Example: nirs.io.saveDotNirs( dod , 'C:\myData' , '%Group%/%index%_Sub%SubjectID%.nirs' );
%   could produce C:\myData\Control\01_Sub1023.nirs, C:\myData\Clinical\02_Sub1401.nirs, etc 
%

%% Initialize and set defaults
if nargin<2, out_dir = ''; end
if nargin<3, filename_pattern = ''; end

if isempty(out_dir), out_dir = pwd; end
out_dir = GetFullPath(out_dir);
if ~exist(out_dir,'dir'), mkdir(out_dir); end

if isempty(filename_pattern)
    filename_pattern = '%index%.nirs';
end

%% Generate index string for filename pattern
num_dec = floor(log10(length(data))) + 1;
index = cell(length(data),1);
for i = 1:length(data)
    index{i} = sprintf(sprintf('%%0%ii',num_dec),i);
end
demo = nirs.createDemographicsTable(data);
demo.index = index;

%% Detect demo variable names in filename
vars = strrep(regexp(filename_pattern,'%\w*%','match'),'%','');
consts = regexp(filename_pattern,'%\w*%','split');

unmatched = setdiff(vars,demo.Properties.VariableNames);
if ~isempty(unmatched)
    error('Unknown variable: %s',unmatched);
end

%% Generate filenames
filenames = cell(length(data),1);
for i = 1:length(data)
    filenames{i} = [out_dir filesep consts{1}];
    for j = 1:length(vars)
        filenames{i} = [filenames{i} demo.(vars{j}){i} consts{j+1}];
    end
end

for i = 1:length(data)
    %% Skip overwrite existing files
    if exist(filenames{i},'file'), warning('File exists: %s. Skipping',filenames{i}); continue; end
    if ~exist(fileparts(filenames{i}),'dir')
        mkdir(fileparts(filenames{i}));
    end
    
    %% Sort probe link
    [~,sorted_inds] = sortrows(data(i).probe.link,[3 1 2]);
    data(i).probe.link = data(i).probe.link(sorted_inds,:);
    data(i).data = data(i).data(:,sorted_inds);
    
    %% Create output .nirs data structure
    odata = struct;
    odata.SD = nirs.util.probe2sd( data(i).probe );
    odata.ml = odata.SD.MeasList;
    odata.d = data(i).data;
    odata.t = data(i).time;
    
    odata.s = false(size(odata.t));
    for c = 1:length(data(i).stimulus.keys)
        stim = data(i).stimulus( data(i).stimulus.keys{c} );
        odata.s = odata.s | stim.getStimVector( odata.t );
    end
    
    %% Write
    save(filenames{i},'-struct','odata');
    
end

end