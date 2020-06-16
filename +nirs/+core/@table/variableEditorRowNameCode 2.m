function [metadataCode,warnmsg] = variableEditorRowNameCode(this,varName,index,rowName)
%   This function is undocumented and will change in a future release

% Generate MATLAB command to modify row names in the row positions 
% defined by the index input.

%   Copyright 2011-2012 The MathWorks, Inc.

warnmsg = '';

% Validation
if isempty(rowName)
     error(message('MATLAB:table:InvalidRowNames'));
end
if ~isempty(this.Properties.RowNames) && any(strcmp(this.Properties.RowNames,rowName))
    error(message('MATLAB:table:DuplicateRowNames',rowName));
end

metadataCode = [varName '.Properties.RowNames{' num2str(index) '} = ''' rowName ''';'];



