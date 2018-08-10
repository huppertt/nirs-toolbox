function [names,indices,classes,iscellStr,charArrayWidths] = variableEditorColumnNames(a)
%   This function is undocumented and will change in a future release

% Undocumented method used by the Variable Editor to determine the names and
% properties of table columns.

%   Copyright 2011-2014 The MathWorks, Inc.

names = a.Properties.VariableNames;

% indices identifies the column positions of each table variable with
% an additional last value of indices is the column after the last column of the
% dataset.
if nargout >= 2
    if isempty(a)
        indices = 1:size(a,2)+1;
    else
        indices = cumsum([1 cellfun(@(x) size(x,2)*ismatrix(x)*~ischar(x)*~isa(x,'dataset')*~isa(x,'table')...
    +ischar(x)+isa(x,'dataset')+isa(x,'table'),a.data)]);
    end
end

if nargout>=3
    classes = cellfun(@class,a.data,'UniformOutput',false); 
end


if nargout>=4
    iscellStr = false(length(names),1);
    a_data = a.data;
    for index=1:length(names)
        if strcmp(classes{index},'cell')
            iscellStr(index) = iscellstr(a_data{index});
        end
    end
end

% Determine the number of chars in any char array columns
if nargout>=5
    charArrayWidths = zeros(length(names),1);
    for index=1:length(names)
        if strcmp(classes{index},'char')
            charArrayWidths(index) = size(a_data{index},2);
        end
    end
end
