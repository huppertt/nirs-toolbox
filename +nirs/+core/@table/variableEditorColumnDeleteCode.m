function [out,warnmsg] = variableEditorColumnDeleteCode(this,varName,colIntervals)
%   This function is undocumented and will change in a future release

% Generate MATLAB command to delete variables in column positions defined
% by the 2-column colIntervals matrix. It is assumed that column intervals
% are disjoint, in monotonic order, and bounded by the number of columns 
% in the dataset array.

%   Copyright 2011-2012 The MathWorks, Inc.

warnmsg = '';
[~,varIndices] = variableEditorColumnNames(this);
varNames = this.Properties.VariableNames;

varNameArray = '';
numVariables = 0;
for k=1:size(colIntervals,1)
    index1 = find(varIndices(1:end-1)<=colIntervals(k,1),1,'last');
    index2 = find(varIndices(1:end-1)<=colIntervals(k,2),1,'last');
    
    % Construct a comma separated list of table variable names
    for index=index1:index2
        varNameArray = [varNameArray '''' varNames{index} '''']; %#ok<AGROW>
        if index<index2 || k<size(colIntervals,1)
            varNameArray = [varNameArray ',']; %#ok<AGROW>
        end
    end
    numVariables = numVariables+index2-index1+1;
end
if numVariables>1
    varNameArray = ['{' varNameArray '}'];
end
out = [varName '(:,' varNameArray ') = [];'];
