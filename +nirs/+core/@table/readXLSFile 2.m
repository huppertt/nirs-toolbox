function t = readXLSFile(xlsfile,args)
%READXLSFILE Read in an XLS file and create a table.

%   Copyright 2012 The MathWorks, Inc.

import matlab.internal.tableUtils.validateLogical

pnames = {'ReadVariableNames' 'ReadRowNames' 'TreatAsEmpty' 'Sheet' 'Range' 'Basic'};
dflts =  {               true          false             {}      ''      ''  false };
[readVarNames,readRowNames,treatAsEmpty,sheet,range,basic] ...
                   = matlab.internal.table.parseArgs(pnames, dflts, args{:});
readRowNames = validateLogical(readRowNames,'ReadRowNames');
readVarNames = validateLogical(readVarNames,'ReadVariableNames');
basic = validateLogical(basic,'Basic');

if isempty(treatAsEmpty)
    treatAsEmpty = cell(0,1);
elseif ischar(treatAsEmpty) && ~isrow(treatAsEmpty)
    % textscan does something a little obscure when treatAsEmpty is char but
    % not a row vector, disallow that here.
    error(message('MATLAB:readtable:InvalidTreatAsEmpty'));
elseif ischar(treatAsEmpty) || iscellstr(treatAsEmpty)
    if ischar(treatAsEmpty), treatAsEmpty = cellstr(treatAsEmpty); end
    % Trim insignificant whitespace to be consistent with what's done for text files.
    treatAsEmpty = strtrim(treatAsEmpty);
    if any(~isnan(str2double(treatAsEmpty))) || any(strcmpi('nan',treatAsEmpty))
        error(message('MATLAB:readtable:NumericTreatAsEmpty'));
    end
else
    error(message('MATLAB:readtable:InvalidTreatAsEmpty'));
end

if basic
    warnState = warning('off','MATLAB:xlsread:Mode');
    cleanup = onCleanup(@() warning(warnState));
    [numeric,txt,raw] = xlsread(xlsfile,sheet,range,'basic');
else
    [numeric,txt,raw] = xlsread(xlsfile,sheet,range);
end
if isempty(txt) && isempty(numeric)
    t = table;
    return
end
clear numeric txt

if readVarNames
    varNames = raw(1,:);
    if ~iscellstr(varNames)
        varNames = convertCell2Str(varNames);
    end
    raw(1,:) = [];
else
    % If reading row names, number remaining columns beginning from 1, we'll drop Var0 below.
    varNames = matlab.internal.table.dfltVarNames((1:size(raw,2))-readRowNames);
end

if readRowNames
    rowNames = raw(:,1);
    if ~iscellstr(rowNames)
        rowNames = convertCell2Str(rowNames);
    end
    dimNames = matlab.internal.table.dfltDimNames;
    if readVarNames, dimNames{1} = varNames{1}; end
    varNames(1) = [];
    raw(:,1) = [];
end

nvars = length(varNames);
t_data = cell(1,nvars);

treatAsEmpty{end+1} = ''; % always treat the empty string as, well, empty
for j = 1:nvars
    var_j = raw(:,j);
    num_j = cellfun(@isnumeric,var_j);
    if all(num_j)
        % Convert columns that are all double to a double var
        t_data{j} = cat(1,var_j{:});
    else
        % Any TreatAsEmpty strings, or any strings of the form 'NaN', in the raw
        % data should be treated as numeric.  TreatAsEmpty has been strtrim'ed,
        % do the same for the data.
        char_j = cellfun('isclass',var_j,'char');
        s = strtrim(var_j(char_j));
        tf = strcmpi('NaN',s);
        for i = 1:length(treatAsEmpty), tf = tf | strcmp(treatAsEmpty{i},s); end
        num_j(char_j) = tf;
        if all(num_j)
            % Convert columns that are all double to a double var
            var_j(num_j & char_j) = {NaN}; % one of the treatAsEmpty strings
            t_data{j} = cat(1,var_j{:});
        elseif all(char_j)
            % Convert columns that are all char to a string var
            t_data{j} = var_j;
        else
            log_j = cellfun(@islogical,var_j);
            if all(num_j | log_j)
                % Convert columns that all logical to a logical var, mixed double/logical to a double var
                var_j(num_j & char_j) = {NaN}; % one of the treatAsEmpty strings
                t_data{j} = cat(1,var_j{:});
            else
                % Convert columns that are mixed anything else to a string var
                t_data{j} = convertCell2Str(var_j);
            end
        end
    end
end

t = table(t_data{:});

% Set the var names.  These will be modified to make them valid, and the
% original strings saved in the VariableDescriptions property.  Fix up
% duplicate or empty names.
t = setVarNames(t,varNames(:)',[],true,true,true);

if readRowNames
    t = setRowNames(t,rowNames,[],true,true); % Fix up duplicate or empty names
    t = setDimNames(t,dimNames,true,true); % Fix up duplicate or empty names
end


%-----------------------------------------------------------------------------
function x = convertCell2Str(x)

for i = 1:length(x)
    y = x{i};
    if isnumeric(y)
        if isnan(y)
            % A numeric NaN means the cell was empty, make that the empty string
            s = '';
        else
            s = num2str(y);
        end
    elseif islogical(y)
        if y
            s = 'true';
        else
            s = 'false';
        end
    elseif ischar(y)
        s = y;
    else
        error(message('MATLAB:readtable:UnexpectedClass', class( y )));
    end
    x{i} = s;
end
