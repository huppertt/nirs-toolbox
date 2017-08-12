function [metadataCode,warnmsg] = variableEditorRowNameCode(this,varName,index,obsName)
%   This function is undocumented and will change in a future release

% Generate MATLAB command to modify observations names in the row positions 
% defined by the index input.

%   Copyright 2012 The MathWorks, Inc.

warnmsg = '';

% Validation
if isempty(obsName)
     error(message('stats:dataset:setobsnames:InvalidObsnames'));
end
if ~isempty(this.Properties.ObsNames) && any(strcmp(this.Properties.ObsNames,obsName))
    error(message('stats:dataset:setobsnames:DuplicateObsnames',obsName));
end

metadataCode = [varName '.Properties.ObsNames{' num2str(index) '} = ''' obsName ''';'];



