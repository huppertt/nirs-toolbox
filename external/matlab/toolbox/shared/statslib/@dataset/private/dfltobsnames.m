function names = dfltobsnames(obsIndices,oneName)
prefix = getString(message('stats:dataset:uistrings:DfltObsNamePrefix'));
if nargin < 2 || ~oneName
    if isempty(obsIndices)
        names = cell(0,1);
    else
        names = strcat({prefix},num2str(obsIndices(:),'%-d'));
    end
else
    names = strcat(prefix,num2str(obsIndices,'%-d'));
end
