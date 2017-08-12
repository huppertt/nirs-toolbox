function [ungroupCode,msg] = variableEditorUngroupCode(this,varName,col)
%   This function is undocumented and will change in a future release

% Generate MATLAB command to ungroup variables in column position defined
% by the input variabe col.

%   Copyright 2011 The MathWorks, Inc.

msg = '';
[~,varIndices] = variableEditorColumnNames(this);
index = find(varIndices(1:end-1)<=col,1,'last');

% If the grouped index is not the first, begin with a subsref for the 
% first index-1 dataset variables.
if index>1 
    ungroupCode = [varName ' = [' varName '(:,1:' num2str(index-1) '),'];
else
    ungroupCode = [varName ' = ['];
end

% Add code to insert a separate dataset for each ungrouped dataset
% variable, e.g., ... dataset(x.GroupHeader(:,1),'VarName','GroupHeader1'),
% ...
datasetVarName = this.Properties.VarNames{index};
if ~isempty(this.Properties.Units)
    datasetVarUnits = this.Properties.Units{index};
else
    datasetVarUnits = '';
end
if ~isempty(this.Properties.VarDescription) 
    datasetVarDescription = this.Properties.VarDescription{index};
else
    datasetVarDescription = '';
end
for k=varIndices(index):varIndices(index+1)-1
    newDatasetVarName = [datasetVarName num2str(num2str(k-varIndices(index)+1))];
    newDatasetVarName = localGetUniqueVarName(newDatasetVarName,this.Properties.VarNames);
    datasetExpression = ['dataset(' varName '.' datasetVarName '(:,' num2str(k-varIndices(index)+1) ...
        '),''VarNames'',''' newDatasetVarName ''')'];
    ungroupCode = [ungroupCode datasetExpression]; %#ok<AGROW>
    if k<varIndices(index+1)-1
        ungroupCode = [ungroupCode ',']; %#ok<AGROW>
    end
end

% If the grouped index is not the last, end with a subsref for the 
% last index+1 dataset variables.
if index<size(this,2)
    ungroupCode = [ungroupCode ',' varName '(:,' num2str(index+1) ':end)];'];
else
    ungroupCode = [ungroupCode '];'];
end

% Copy any non-empty metadata onto the ungrouped dataset variables.
newIndices = varIndices(index+1)-varIndices(index);
if ~isempty(datasetVarUnits)
    for k=index:index+newIndices-1
        ungroupCode = [ungroupCode varName '.Properties.Units{' num2str(k) '} = ''' datasetVarUnits ''';']; %#ok<AGROW>
    end
end
if ~isempty(datasetVarDescription)
    for k=index:index+newIndices-1
        ungroupCode = [ungroupCode varName '.Properties.VarDescription{' num2str(k) '} = ''' datasetVarDescription ''';']; %#ok<AGROW>
    end
end

function varName = localGetUniqueVarName(newVarName,varNames)

% Append digits to guarantee a unique variable name
ind = 1;
varName = newVarName;
while any(strcmp(varName,varNames))
    varName = [newVarName,num2str(ind)];
    ind=ind+1;
end

%hospital = [hospital(:,1:5),dataset(hospital.BloodPressure(:,1),'VarNames',{'BloodPressure_1'}),dataset(hospital.BloodPressure(:,2),'VarNames',{'BloodPressure_2'}),hospital(:,7:end)];"