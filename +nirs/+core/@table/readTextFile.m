function t = readTextFile(file, args)
%READFILE Read in a delimited text file and create a table.

%   Copyright 2012-2014 The MathWorks, Inc.
import matlab.internal.tableUtils.validateLogical

pnames = {'ReadVariableNames' 'ReadRowNames' 'Delimiter' 'Format' 'TreatAsEmpty' 'HeaderLines' 'FileEncoding'};
dflts =  {               true          false     'comma'       []             {}            0          ''};
[readVarNames,readRowNames,delimiter,format,treatAsEmpty,headerLines, fileEncoding, supplied,otherArgs] ...
    = matlab.internal.table.parseArgs(pnames, dflts, args{:});

readRowNames = validateLogical(readRowNames,'ReadRowNames');
readVarNames = validateLogical(readVarNames,'ReadVariableNames');

inferFormat = ~supplied.Format;

tab = sprintf('\t');

% Set the delimiter, recognizing aliases
if ischar(delimiter)
    whiteSpace = sprintf(' \b\t');
    switch delimiter
        case {'tab', '\t', tab}
            delimiter = tab;
            whiteSpace = sprintf(' \b');
        case {'space',' '}
            delimiter = ' ';
            whiteSpace = sprintf('\b\t');
        case {'comma', ','}
            delimiter = ',';
        case {'semi', ';'}
            delimiter = ';';
        case {'bar', '|'}
            delimiter = '|';
        otherwise
            if ~isscalar(delimiter)
                error(message('MATLAB:readtable:InvalidDelimiter'));
            end
    end
else
    error(message('MATLAB:readtable:InvalidDelimiter'));
end

if isempty(treatAsEmpty)
    treatAsEmpty = {};
elseif ischar(treatAsEmpty) && ~isrow(treatAsEmpty)
    % textscan does something a little obscure when treatAsEmpty is char but
    % not a row vector, disallow that here.
    error(message('MATLAB:readtable:InvalidTreatAsEmpty'));
elseif ischar(treatAsEmpty) || iscellstr(treatAsEmpty)
    % TreatAsEmpty is only ever applied to numeric fields in the file, and
    % textscan ignores leading/trailing whitespace for those fields, so trim
    % insignificant whitespace.
    treatAsEmpty = strtrim(treatAsEmpty);
    if any(~isnan(str2double(treatAsEmpty))) || any(strcmpi('nan',treatAsEmpty))
        error(message('MATLAB:readtable:NumericTreatAsEmpty'));
    end
else
    error(message('MATLAB:readtable:InvalidTreatAsEmpty'));
end

if ~isscalar(headerLines) || ~isnumeric(headerLines) || ...
        (headerLines < 0) || (round(headerLines) ~= headerLines)
    error(message('MATLAB:readtable:InvalidHeaderLines'));
end

% textscan has many parameters that the local fread-based function does not
% support. If a format is not specified, may end up falling back to that code,
% so can't accept those params unless a forrmat is provided.
if inferFormat && ~isempty(otherArgs)
    try
        % Give error based on whether we have a textscan param, or something
        % completely unknown.
        textscan('a','%s',otherArgs{:});
    catch ME
        if strcmp(ME.identifier,'MATLAB:textscan:UnknownOption')
            error(message('MATLAB:table:parseArgs:BadParamName',otherArgs{1}));
        end
    end
    error(message('MATLAB:readtable:MustProvideFormat',otherArgs{1}));
end

% Open the file.
fid = fopen(file,'rt','n',fileEncoding); % text mode: CRLF -> LF on windows (no-op on linux)
if fid == -1
    % Try again with default extension if there wasn't one
    [~,~,fx] = fileparts(file);
    if isempty(fx)
        fid = fopen([file '.txt'],'rt','n',fileEncoding); % text mode: CRLF -> LF on windows (no-op on linux)
    end
end
if fid == -1
    error(message('MATLAB:readtable:OpenFailed',file));
end
cleanFid = onCleanup(@()fclose(fid));

if readVarNames
    % Read in the first line of var names as a single string, skipping any leading
    % blank lines.
    raw = textscan(fid, '%[^\r\n]', 1, 'Whitespace', '', 'headerlines', headerLines, otherArgs{:});
    if isempty(raw{1}) || isempty(raw{1}{1})
        if inferFormat
            t = table.empty(0,0);
            return
        else
            vnline = ' '; % whitespace
        end
    else
        vnline = raw{1}{1};
    end
else
    % Skip header lines
    textscan(fid, '%[^\r\n]', 0, 'Whitespace', '', 'headerlines', headerLines, otherArgs{:});
    
    vnline = '';
end

% Guess at a format string for the data.
numNonDataLines = headerLines + readVarNames;
if inferFormat
    
    % Save current position in file
    startDataPosn = ftell(fid);
    
    % Guess at format string based on first line of data
    format = guessformat(fid, vnline, delimiter, whiteSpace, treatAsEmpty);
    
    % Rewind fid to saved position
    fseek(fid, startDataPosn, 'bof');
    
    % Read in the data using the guessed-at format.  textscan
    % automatically skips blank lines.  Even if there's nothing left in the file,
    % textscan will return the right types in raw.
    raw = textscan(fid, format, ...
        'Delimiter', delimiter, 'Whitespace', whiteSpace, ...
        'TreatAsEmpty', treatAsEmpty, 'EndOfLine', '\r\n');
    
    if ~feof(fid)
        % textscan failed because some column had the "wrong" value in it, i.e.,
        % not the same type as in that column in the first data line.
        
        % Rewind fid to saved position
        fseek(fid, startDataPosn, 'bof');
        
        % Try guessing the format by going through all lines of the file
        format = updateformat(fid, format, delimiter, whiteSpace, treatAsEmpty, numNonDataLines);
        
        raw = textscan(fid, format, ...
            'Delimiter', delimiter, 'Whitespace', whiteSpace, ...
            'TreatAsEmpty', treatAsEmpty, 'EndOfLine', '\r\n');
        
        if ~feof(fid)
            % The format found by updatewholeformat is not correct.
            % This will only happen in very specific cases, here is
            % an example:
            %    1, 11, aa, 21, xx
            %    2, 12  bb, 22, yy,
            %    3, EE, cc, 23, zz
            % The problem occurs because the comma between 12 and bb is
            % missing *and* a comma is added at the end of the line.
            
            % We set the format of all variables to %q.
            format(2:2:end) = 'q';
            
            fseek(fid,startDataPosn,'bof');
            raw = textscan(fid, format, ...
                'Delimiter', delimiter, 'Whitespace', whiteSpace, ...
                'TreatAsEmpty', treatAsEmpty, 'EndOfLine', '\r\n');
        end
    end
    
else % ~inferFormat
    
    % Read in the data using the specified format.  textscan
    % automatically skips blank lines.  Even if there is nothing left in the file,
    % textscan will return the right types in raw.
    raw = textscan(fid, format, ...
        'Delimiter', delimiter, 'Whitespace', whiteSpace, 'TreatAsEmpty', treatAsEmpty, otherArgs{:});
    
    if ~feof(fid)
        m = message('MATLAB:readtable:CouldNotReadEntireFileWithFormat');
        baseME = MException(m.Identifier,'%s',getString(m));
        % If all the cells in raw are the same length, textscan stopped at the
        % start of a line.  If the first few cells have length L, and the rest
        % L-1, then textscan stopped mid-line.  We can be helpful here.
        varlens = cellfun(@(x)size(x,1),raw);
        dvarlens = diff(varlens);
        locs = find(dvarlens);
        if isempty(locs) || (isscalar(locs) && dvarlens(locs)==-1)
            errLine = min(varlens) + 1 + numNonDataLines;
            m = message('MATLAB:readtable:ReadErrorOnLine',errLine);
            throw(addCause(baseME,MException(m.Identifier,'%s',getString(m))));
            % Otherwise, something else happened for which we have no specific advice.
        else
            throw(baseME);
        end
    end
end

if readVarNames
    % determine variable names from the variable names line
    varNames = ...
        matlab.internal.table.determineVarNames(vnline,format,delimiter,whiteSpace,true,otherArgs);
    
    if numel(varNames) ~= length(raw)
        error(message('MATLAB:readtable:ReadVarNamesFailed',file,length(raw),numel(varNames)));
    end
else
    % If reading row names, number remaining columns beginning from 1, we'll drop Var0 below.
    varNames = matlab.internal.table.dfltVarNames((1:length(raw))-readRowNames);
end

if isempty(raw) % i.e., if the file had no data
    t_data = cell(length(raw));
else
    varlen = unique(cellfun(@(x)size(x,1),raw));
    if ~isscalar(varlen)
        if inferFormat
            error(message('MATLAB:readtable:UnequalVarLengthsFromFileNoFormat'));
        else
            error(message('MATLAB:readtable:UnequalVarLengthsFromFileWithFormat'));
        end
    end
    
    if readRowNames
        rowNames = raw{1};
        if ischar(rowNames)
            rowNames = cellstr(rowNames);
        elseif isnumeric(rowNames)
            rowNames = cellstr(num2str(rowNames));
        elseif ~iscellstr(rowNames)
            error(message('MATLAB:readtable:RowNamesVarNotString', class(rowNames)));
        end
        raw(1) = [];
        dimNames = matlab.internal.table.dfltDimNames;
        if readVarNames, dimNames{1} = varNames{1}; end
        varNames(1) = [];
    end
    t_data = raw(:)';
end

t = table(t_data{:});

% Set the var names.  These will be modified to make them valid, and the
% original strings saved in the VariableDescriptions property.  Fix up
% duplicate or empty names.
t = setVarNames(t,varNames(:)',[],true,true,true);

if ~isempty(raw) && readRowNames
    t = setRowNames(t,rowNames,[],true,true); % Fix up duplicate or empty names
    t = setDimNames(t,dimNames,true,true); % Fix up duplicate or empty names
end

% Guess at a format string for the data, based on the first line of data
% accessible to fid and the variable names (if vnline isn't empty)
function format = guessformat(fid, vnline, delimiter, whiteSpace, treatAsEmpty)
% Read in the first line of data as a single string, skipping any leading
% blank lines.  Then back up.
startDataPosn = ftell(fid);
raw = textscan(fid, '%[^\r\n]', 1, 'Whitespace', whiteSpace);
fseek(fid,startDataPosn,'bof');

% If first line is empty, use vnline (variable names) instead
if isempty(raw{1}) || isempty(raw{1}{1})
    if ~isempty(vnline)
        nvars = nnz(vnline==delimiter) + 1;
    else
        nvars = 0;
    end
    if nvars > 0
        format = repmat('%f',1,nvars);
    else
        format = '%*s'; % textscan does not accept an empty format
    end
else
    % Determine how many data columns there are.
    format = matlab.internal.table.determineFormatString(raw{1}{1}, ...
        delimiter, whiteSpace, treatAsEmpty);    
end

% Guess at a format string for the data, based on reading all lines
% of the file using textscan, and on the original format string
% from the preceding call to guessformat.
function format = updateformat(fid, format, delimiter, whiteSpace, treatAsEmpty, numNonDataLines)

nvars = numel(format)/2;

% format_parse: textscan reads entries, but generates no output.
% This is a faster way to determine if there is a problem on a given line.
format_parse = repmat('%*f', 1, nvars);
format_parse(3:3:end) = format(2:2:end);

% Save current position in file
startDataPosn = ftell(fid);

% Number of lines to be read with textscan at once
block_size = 1e2;

while(~feof(fid))
    
    % Read block_size lines from file as string
    rawline = textscan(fid, '%[^\r\n]', block_size, 'Whitespace', '', 'EndOfLine', '\r\n');
    
    if isempty(rawline{1}) % reached end of file
        break;
    end
    rawline = rawline{1};
    
    % Go through each line
    for jj=1:numel(rawline)
        
        % Test that file is rectangular (nbr. of delimiters per line
        % doesn't change)
        numDelimiter = nnz(rawline{jj}==delimiter);
        if numDelimiter ~= nvars-1
            errLine = jj + numNonDataLines; % account for nonDataLines not read into rawline
            error(message('MATLAB:readtable:BadFileFormat', errLine, errLine, numDelimiter, nvars-1));
        end
        
        % Try parsing line
        [~,pos] = textscan(rawline{jj}, format_parse, 1, ...
            'Delimiter', delimiter, 'Whitespace', whiteSpace, 'TreatAsEmpty', treatAsEmpty);
        
        % If parsing line failed, update the format string. We pass an
        % empty cell as an argument to determineFormatString as readtable
        % currently does not allow passing additional textscan arguments
        % for format detection
        if pos ~= length(rawline{jj})
            format = matlab.internal.table.determineFormatString(rawline{jj}, ...
                delimiter, whiteSpace, treatAsEmpty, {}, format);
            format_parse(3:3:end) = format(2:2:end);
        end
    end
end

% Rewind fid to saved position
fseek(fid,startDataPosn,'bof');

