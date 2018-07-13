function currFormat = variableEditorMetadata(this)
% This function is undocumented and will change in a future release

% Retrieves the format for any datetime columns in the table, which
% is needed for the variable editor.

% Copyright 2014 The MathWorks, Inc.

% Get the Format for any datetime columns
[~, varIndices, colClasses] = variableEditorColumnNames(this);
idx = strcmp(colClasses, 'datetime');
currFormat = '';
for col=varIndices(idx)
    % For table we currently don't provide a mechanism to display/choose
    % multiple formats, so we will just return the first one found.
    d = this.data{col};
    if ~strcmp(d.TimeZone, 'UTCLeapSeconds')
        currFormat = d.Format;
        break;
    end
end

