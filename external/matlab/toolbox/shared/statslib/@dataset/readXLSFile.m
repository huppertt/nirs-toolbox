function a = readXLSFile(a,xlsfile,args)
%READXLSFILE Read in an XLS file and create a dataset array.

%   Copyright 2011-2012 The MathWorks, Inc.


pnames = {'readvarnames' 'readobsnames' 'sheet' 'range'};
dflts =  {          true          false      ''      ''};
[readvarnames,readobsnames,sheet,range] ...
                   = dataset.parseArgs(pnames, dflts, args{:});
readobsnames = onOff2Logical(readobsnames,'ReadObsNames');
readvarnames = onOff2Logical(readvarnames,'ReadVarNames');

[numeric,txt,raw] = xlsread(xlsfile,sheet,range);
if isempty(numeric) && isempty(txt)
    return
end
clear numeric txt

if readvarnames
    varnames = raw(1,:);
    if ~iscellstr(varnames)
        varnames = cellfun(@convertCell2Str, varnames, 'UniformOutput',false);
    end
    raw(1,:) = [];
else
    varnames = dfltvarnames(1:size(raw,2));
end

a.nobs = size(raw,1);
if readobsnames
    obsnames = raw(:,1);
    if ~iscellstr(obsnames)
        obsnames = cellfun(@convertCell2Str, obsnames, 'UniformOutput',false);
    end
    dimnames = a.props.DimNames;
    dimnames{1} = varnames{1};
    varnames(1) = [];
    raw(:,1) = [];
    a = setobsnames(a,obsnames);
    a = setdimnames(a,dimnames);
end

a.nvars = length(varnames);
a.data = cell(1,a.nvars);

% Set the var names.  These will be modified to make them valid, and the
% original strings saved in the VarDescription property.  Allow empty and
% duplicate names.
a = setvarnames(a,varnames(:)',[],true,true,true);

for j = 1:a.nvars
    var_j = raw(:,j);
    num_j = cellfun(@isnumeric,var_j);
    if all(num_j)
        % convert columns that are all double to a double var
        a.data{j} = cell2mat(var_j);
    else
        char_j = cellfun(@ischar,var_j);
        if all(char_j)
            % convert columns that are all char to a string var
            a.data{j} = var_j;
        else
            log_j = cellfun(@islogical,var_j);
            if all(log_j)
                % convert columns that are all logical to a logical var
                a.data{j} = cell2mat(var_j);
            elseif all(num_j | log_j)
                % convert columns that are mixed double/logical to a double var
                a.data{j} = cellfun(@double,var_j);
            else
                % convert columns that are mixed anything else to a string var
                a.data{j} = cellfun(@convertCell2Str, var_j, 'UniformOutput',false);
            end
        end
    end
end


%-----------------------------------------------------------------------------
function s = convertCell2Str(n)
if isnumeric(n)
    if isnan(n)
        % A numeric NaN means the cell was empty, make that the empty string
        s = '';
    else
        s = num2str(n);
    end
elseif islogical(n)
    if n
        s = 'true';
    else
        s = 'false';
    end
elseif ischar(n)
    s = n;
else
    error(message('stats:dataset:dataset:UnexpectedClass', class( n )));
end
    
