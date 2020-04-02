function [tbdata, run] = readEPrimePDAT(file, path)
% 
% function tbdata = readEPrimePDAT(file, path)
% 
% Returns table data for getStimDesign GUI from an E-Prime PDAT file.
% 
% Davneet Minhas
% June 4, 2010
% 

% Set Number of Columns in PDAT File
n = '';
for i=1:21
    n = [n,'%s'];
end

% Read File into Cell String
fid = fopen([path, '/', file]);
data = textscan(fid, n, 'delimiter', '\t', 'CollectOutput', 1);
all = data{1};
fclose(fid);

% Check for non-Strings (Very Inefficient Implementation)
[row, col] = size(all);
for r = 1:row
    for c = 1:col
        if ~isstr(all{r,c})
            all{r,c} = '';
        end
    end
end

% Determine Relevant Columns
[rRun, cRun] = find(ismember(all, 'Run'));
[rOns, cOns] = find(ismember(all, 'OnsetTime'));
[rDur, cDur] = find(ismember(all, 'Duration'));
[rCon, cCon] = find(ismember(all, 'ConditionId'));

% Determine Number of Runs
pr = unique(all(rRun+1:end,cRun));

% Select Run Number
run = inputdlg(['Select Run # (', pr{1}, '-', pr{end}, ')']);

% Isolate Run Number Data
col = find(ismember(all(:,cRun), run{1}));
col = col(col > rRun);
rundata = all(col,:);

% Create Table
tbdata = cell(size(rundata, 1), 4);
zero = str2num(rundata{1,cOns})/1000;
for i=1:size(rundata,1)
    tbdata{i,1} = str2num(rundata{i,cOns})/1000;
    tbdata{i,2} = str2num(rundata{i,cDur})/1000;
    tbdata{i,3} = 1;
    tbdata{i,4} = rundata{i,cCon};
    
    % Set Initial Onset to Zero    
    tbdata{i,1} = tbdata{i,1} - zero;
end


end