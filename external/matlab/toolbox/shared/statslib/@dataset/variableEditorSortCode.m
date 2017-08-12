function [sortCode,msg] = variableEditorSortCode(~,varName,datasetVariableNames,direction)
%   This function is undocumented and will change in a future release

% Generate MATLAB command to sort dataset observations. The direction input
% is true for ascending sorts, false otherwise.

%   Copyright 2011 The MathWorks, Inc.

msg = '';
if iscell(datasetVariableNames)
    datasetVariableNameString = '{';
    for k=1:length(datasetVariableNames)-1
        datasetVariableNameString = [datasetVariableNameString datasetVariableNames{k} ',']; %#ok<AGROW>
    end
    datasetVariableNamesString = [datasetVariableNameString datasetVariableNames{end} '}'];
else
    datasetVariableNamesString = datasetVariableNames;
end
if direction
    sortCode = [varName ' = sortrows(' varName ',' datasetVariableNamesString ',''ascend'');'];
else
    sortCode = [varName ' = sortrows(' varName ',' datasetVariableNamesString ',''descend'');'];
end



