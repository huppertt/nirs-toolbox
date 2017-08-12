function a = fillInDataset(a,newLen,newObsNames)
% Fill in variables that are too short in a dataset
for j = 1:a.nvars
    if size(a.data{j},1) < newLen
        a.data{j} = lengthenVar(a.data{j}, newLen);
    end
end

% If the original dataset had observation names, append the names for the new
% observations, or append default names
if ~isempty(a.obsnames)
    if ~isempty(newObsNames)
        a.obsnames = [a.obsnames; newObsNames];
    elseif newLen > a.nobs
        a.obsnames = [a.obsnames; dfltobsnames((a.nobs+1):newLen)];
    end
% If the new observations have observation names and the original dataset
% doesn't, create default names for the dataset
elseif ~isempty(newObsNames) % && isempty(a.obsnames)
    if a.nobs > 0
        a.obsnames = [dfltobsnames(1:a.nobs); newObsNames];
    else
        a.obsnames = newObsNames;
    end
end
% Otherwise, do not create observation names if there were none.

a.nobs = newLen;
