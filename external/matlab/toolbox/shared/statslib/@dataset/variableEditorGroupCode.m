function [groupCode,msg] = variableEditorGroupCode(this,varName,startCol,endCol)
%   This function is undocumented and will change in a future release

% Generate MATLAB command to group variables in column range startCol to
% endCol

%   Copyright 2011-2013 The MathWorks, Inc.

msg = '';
[varNames,varIndices] = variableEditorColumnNames(this);
startIndex = find(varIndices(1:end-1)<=startCol,1,'last');
endIndex = find(varIndices(1:end-1)<=endCol,1,'last');

% If startIndex is not the first, begin with a subsref for the 
% first index-1 dataset variables.
if startIndex>1 
    groupCode = [varName ' = [' varName '(:,1:' num2str(startIndex-1) '),dataset(['];
else
    groupCode = [varName ' = [dataset(['];
end

% Add code to insert a dataset for the grouped columns
% e.g., dataset([x.variableName1,x.variableName2,...])
groupedVarName = '';
for k=startIndex:endIndex
    groupedVarName = sprintf('%s%s',groupedVarName,varNames{k});
    groupCode = sprintf('%s%s.%s',groupCode,varName,varNames{k});
    if k<endIndex
       groupCode = sprintf('%s,',groupCode);
       groupedVarName = sprintf('%s_',groupedVarName);
    end   
end
groupCode = [groupCode '], ''varNames'',''' matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(groupedVarName), {}, namelengthmax) ''')'];

% Add the part of the dataset after the grouped dataset variables
if endIndex<size(this,2)
    groupCode = [groupCode ',' varName '(:,' num2str(endIndex+1) ':end)];'];
else
    groupCode = [groupCode '];'];
end