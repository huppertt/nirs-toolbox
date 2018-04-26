function [ungroupCode,msg] = variableEditorUngroupCode(this,varName,col)
%   This function is undocumented and will change in a future release

% Generate MATLAB command to ungroup variables in column position defined
% by the input variabe col.

%   Copyright 2011-2012 The MathWorks, Inc.

msg = '';
[~,varIndices] = variableEditorColumnNames(this);
index = find(varIndices(1:end-1)<=col,1,'last');

% If the grouped index is not the first, begin with a subsref for the 
% first index-1 table variables.
if index>1 
    ungroupCode = [varName ' = [' varName '(:,1:' num2str(index-1) '),'];
else
    ungroupCode = [varName ' = ['];
end

% Add code to insert a separate table for each ungrouped table
% variable, e.g., ... table(x.GroupHeader(:,1),'VariableNames',{'GroupHeader1'}),
% ...
tableVarName = this.Properties.VariableNames{index};
if ~isempty(this.Properties.VariableUnits)
    tableVarUnits = this.Properties.VariableUnits{index};
else
    tableVarUnits = '';
end
if ~isempty(this.Properties.VariableDescriptions) 
    tableVarDescription = this.Properties.VariableDescriptions{index};
else
    tableVarDescription = '';
end
for k=varIndices(index):varIndices(index+1)-1
    newTableVarName = [tableVarName num2str(num2str(k-varIndices(index)+1))];
    newTableVarName = localGetUniqueVarName(newTableVarName,this.Properties.VariableNames);
    tableExpression = ['table(' varName '.' tableVarName '(:,' num2str(k-varIndices(index)+1) ...
        '),''VariableNames'',{''' newTableVarName '''})'];
    ungroupCode = [ungroupCode tableExpression]; %#ok<AGROW>
    if k<varIndices(index+1)-1
        ungroupCode = [ungroupCode ',']; %#ok<AGROW>
    end
end

% If the grouped index is not the last, end with a subsref for the 
% last index+1 table variables.
if index<size(this,2)
    ungroupCode = [ungroupCode ',' varName '(:,' num2str(index+1) ':end)];'];
else
    ungroupCode = [ungroupCode '];'];
end

% Copy any non-empty metadata onto the ungrouped table variables.
newIndices = varIndices(index+1)-varIndices(index);
if ~isempty(tableVarUnits)
    for k=index:index+newIndices-1
        ungroupCode = [ungroupCode varName '.Properties.VariableUnits{' num2str(k) '} = ''' tableVarUnits ''';']; %#ok<AGROW>
    end
end
if ~isempty(tableVarDescription)
    for k=index:index+newIndices-1
        ungroupCode = [ungroupCode varName '.Properties.VariableDescriptions{' num2str(k) '} = ''' tableVarDescription ''';']; %#ok<AGROW>
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

%hospital = [hospital(:,1:5),table(hospital.BloodPressure(:,1),'VariableNames',{'BloodPressure_1'}),table(hospital.BloodPressure(:,2),'VariableNames',{'BloodPressure_2'}),hospital(:,7:end)];"
