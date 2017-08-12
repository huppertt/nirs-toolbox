function [metadataCode,warnmsg] = variableEditorMetadataCode(this,varName,index,propertyName,propertyString)
%   This function is undocumented and will change in a future release

% Generate MATLAB command to modify dataset metadata at positions defined
% by index input. 

%   Copyright 2011-2013 The MathWorks, Inc.

warnmsg = '';
if strcmpi('varnames',propertyName)
    % Validation
    if ~isvarname(propertyString)
         error(message('MATLAB:codetools:InvalidVariableName',propertyString));
    end
    if any(strcmp(this.Properties.VarNames,propertyString))
        error(message('stats:dataset:setvarnames:DuplicateVarnames',propertyString));
    end
    metadataCode = [varName '.Properties.VarNames{' num2str(index) '} = ''' fixquote(propertyString) ''';'];
elseif strcmpi('units',propertyName) || strcmpi('vardescription',propertyName)
    metadataCode = [varName '.Properties.' propertyName '{' num2str(index) '} = ''' fixquote(propertyString) ''';'];
end



