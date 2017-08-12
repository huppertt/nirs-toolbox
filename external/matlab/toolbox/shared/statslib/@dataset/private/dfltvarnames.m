function names = dfltvarnames(varIndices,oneName)
prefix = getString(message('stats:dataset:uistrings:DfltVarNamePrefix'));
if nargin < 2 || ~oneName
    if isempty(varIndices)
        names = cell(1,0);
    else
        names = strcat({prefix},num2str(varIndices(:),'%-d'))';
    end
else
    names = strcat(prefix,num2str(varIndices,'%-d'));
end
