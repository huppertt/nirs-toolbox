function [names,indices,classes,iscellStr,charArrayWidths] = variableEditorColumnNames(a)
%   This function is undocumented and will change in a future release

% Undocumented method used by the Variable Editor to determine the names and
% properties of dataset columns.

%   Copyright 2011-2012 The MathWorks, Inc.

names = a.Properties.VarNames;

% indices identifies the column positions of each dataset variable with
% an additional last value of indices is the column after the last column of the
% dataset. char arrays always occupy one column.
if nargout>=2
    if(isempty(a))
        indices = 1:size(a,2)+1;
    else
        indices = cumsum([1 datasetfun(@(x) size(x,2)*ismatrix(x)*~ischar(x)*~isa(x,'dataset')*~isa(x,'table')...
    +ischar(x)+isa(x,'dataset')+isa(x,'table'),a)]);
    end
end

if nargout>=3
    classes = datasetfun(@class,a,'UniformOutput',false); 
end

% Identify cell array columns
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
