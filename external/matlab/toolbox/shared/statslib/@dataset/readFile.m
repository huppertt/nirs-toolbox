function a = readFile(a,file,args)
%READFILE Read in a delimited text file and create a dataset array.

%   Copyright 2011-2012 The MathWorks, Inc.


pnames = {'readvarnames' 'readobsnames' 'delimiter' 'format' 'treatasempty' 'headerlines'};
dflts =  {          true          false          {}       []             {}            0 };
[readVarNames,readObsNames,delimiter,format,treatAsEmpty,headerLines,supplied,otherArgs] ...
                   = dataset.parseArgs(pnames, dflts, args{:});

readObsNames = onOff2Logical(readObsNames,'ReadObsNames');
readVarNames = onOff2Logical(readVarNames,'ReadVarNames');

inferFormat = ~supplied.format;

tab = sprintf('\t');
lf = sprintf('\n');

% Set the default delimiter, or recognize aliases
if ~supplied.delimiter
    % This is for backwards compatibility.
    if inferFormat
        delimiter = 'tab'; % the tdfread default
    else
        delimiter = ' ';   % the textscan default
    end
end
whiteSpace = sprintf(' \b\t');
switch delimiter
case {'tab', '\t'}
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
    delimiter = delimiter(1);
    warning(message('stats:dataset:dataset:UnrecognizedDelimiter', delimiter(1)));
end

if isempty(treatAsEmpty)
    treatAsEmpty = {};
elseif ischar(treatAsEmpty) && ~isrow(treatAsEmpty)
    % textscan does something a little obscure when treatAsEmpty is char but
    % not a row vector, disallow that here.
    error(message('stats:dataset:dataset:InvalidTreatAsEmpty'));
elseif ischar(treatAsEmpty) || iscellstr(treatAsEmpty)
    treatAsEmpty = strtrim(treatAsEmpty);
    if any(~isnan(str2double(treatAsEmpty))) || any(strcmpi('nan',treatAsEmpty))
        error(message('stats:dataset:dataset:NumericTreatAsEmpty'));
    end
else
    error(message('stats:dataset:dataset:InvalidTreatAsEmpty'));
end

if ~isscalar(headerLines) || ~isnumeric(headerLines) || ...
                (headerLines < 0) || (round(headerLines) ~= headerLines)
    error(message('stats:dataset:dataset:InvalidHeaderLines'));
end
tdfreadHeaderLines = headerLines + readVarNames; % the tdfread local never reads var names

if inferFormat && ~isempty(otherArgs)
    error(message('stats:dataset:parseArgs:BadParamName',otherArgs{1}));
end

% Open the file.
fid = fopen(file,'rt'); % text mode: CRLF -> LF on windows (no-op on linux)
if fid == -1
    error(message('stats:dataset:dataset:OpenFailed',file));
end

if readVarNames
    % Read in the first line of var names as a single string, skipping any leading
    % blank lines.
    raw = textscan(fid,'%s',1,'whitespace',lf,'headerlines',headerLines,otherArgs{:});
    headerLines = 0; % just skipped them
    if isempty(raw{1}) || isempty(raw{1}{1})
        if inferFormat
            fclose(fid);
            a = dataset.empty(0,0);
            return
        else
            vnline = ' '; % whitespace
        end
    else
        vnline = raw{1}{1};
    end
end

% Guess at a format string for the data, based on the first line of data.
if inferFormat
    % Read in the first line of data as a single string, skipping any leading
    % blank lines.  Then back up.
    startDataPosn = ftell(fid);
    raw = textscan(fid,'%s',1,'whitespace',lf,'headerlines',headerLines);
    fseek(fid,startDataPosn,'bof');
    if isempty(raw{1}) || isempty(raw{1}{1})
        if readVarNames
            nvars = length(find(vnline==delimiter)) + 1;
        else
            nvars = 0;
        end
        format = repmat('%f',1,nvars);
    else
        % Determine how many data columns there are.
        line1 = raw{1}{1};
        nvars = length(find(line1==delimiter)) + 1;

        % Read the values in the first line of data as strings, then try to
        % convert each one to numeric.
        strs = textscan(line1,repmat('%q',1,nvars),1,'delimiter',delimiter,'whitespace',whiteSpace);
        format = repmat('%f',1,nvars);
        for i = 1:nvars
            cstr = strs{i};
            if isempty(cstr)
                str = '';
            else
                str = cstr{1};
            end
            num = str2double(str);
            if isnan(num)
                % If the result was NaN, figure out why.
                if isempty(str) || strcmpi(str,'nan') || any(strcmp(str,treatAsEmpty))
                    % NaN came from a nan string, and empty field, or one of the
                    % treatAsEmpty strings, treat this column as numeric.  Note that
                    % because we check str against treatAsEmpty, the latter only works
                    % on numeric columns.  That is what textscan does.
                    format(2*i) = 'f';
                else
                    % NaN must have come from a failed conversion, treat this column
                    % as strings.
                    format(2*i) = 'q';
                end
            else
                % Otherwise the conversion succeeded, treat this column as numeric.
                format(2*i) = 'f';
            end
        end
    end
end

% Read in the data using the guessed-at or specified format.  textscan
% automatically skips blank lines.  Even if there's nothing left in the file,
% textscan will return the right types in raw.
raw = textscan(fid,format,'delimiter',delimiter,'whitespace',whiteSpace, ...
    'headerlines',headerLines,'treatasempty',treatAsEmpty,otherArgs{:});
atEOF = feof(fid);
fclose(fid);

if ~atEOF
    % textscan failed because some column had the "wrong" value in it, i.e.,
    % not the same type as in that column in the first data line.  If no
    % format was given, fall back to the slower but more forgiving TDFREAD.
    % Otherwise, give up.
    if inferFormat
        raw = tdfread(file,delimiter,tdfreadHeaderLines,treatAsEmpty);
    else
        error(message('stats:dataset:dataset:UnequalVarLengthsFromFileWithFormat'));
    end
end

if readVarNames
    % Search for any of the allowable conversions: a '%', followed
    % optionally by '*', followed optionally by 'nnn' or by 'nnn.nnn',
    % followed by one of the type specifiers or a character list or a
    % negative character list.  Keep the '%' and the '*' if it's there,
    % but replace everything else with the 'q' type specifier.
    specifiers = '(n|d8|d16|d32|d64|d|u8|u16|u32|u64|u|f32|f64|f|s|q|c|\[\^?[^\[\%]*\])';
    vnformat = regexprep(format,['\%([*]?)([0-9]+(.[0-9]+)?)?' specifiers],'%$1q');
    % Append a delimiter in case the last varname is missing, so we get a cell
    % containing an empty string, not an empty cell.
    vnline(end+1) = delimiter;
    varNames = textscan(vnline,vnformat,1,'delimiter',delimiter,'whitespace',whiteSpace,otherArgs{:});
    % If a textscan was unable to read a varname, the corresponding cell
    % contains an empty cell.  Remove those.  This happens when trying to
    % put delimiters in the format string instead of using the 'Delimiter'
    % input, because they're just read as part of the string.  This is
    % different than if there is no var name -- in that case, the cell will
    % contain the empty string and be left alone.  setvarnames fixes those
    % later on.
    varNames = varNames(~cellfun('isempty',varNames));
    % Each cell in varnames contains another 1x1 cell containing a string,
    % get those out.  textscan ignores leading whitespace, remove trailing
    % whitspace here.
    varNames = deblank(cellfun(@(c) c{1}, varNames,'UniformOutput',false));
    
    if numel(varNames) ~= length(raw)
        error(message('stats:dataset:dataset:ReadVarnamesFailed',file));
    end
else
    varNames = dfltvarnames(1:length(raw));
end
    
if isempty(raw) % i.e., if the file was empty
    a.nobs = 0;
else
    varlen = unique(cellfun(@(x)size(x,1),raw));
    if ~isscalar(varlen)
        if inferFormat
            error(message('stats:dataset:dataset:UnequalVarLengthsFromFileNoFormat'));
        else
            error(message('stats:dataset:dataset:UnequalVarLengthsFromFileWithFormat'));
        end
    end
    
    a.nobs = varlen;
    if readObsNames
        obsNames = raw{1};
        if ischar(obsNames)
            obsNames = cellstr(obsNames);
        elseif isnumeric(obsNames)
            obsNames = cellstr(num2str(obsNames));
        elseif ~iscellstr(obsNames)
            error(message('stats:dataset:dataset:ObsnamesVarNotString', class(obsNames)));
        end
        raw(1) = [];
        dimnames = a.props.DimNames;
        dimnames{1} = varNames{1};
        varNames(1) = [];
        a = setobsnames(a,obsNames);
        a = setdimnames(a,dimnames);
    end
end
a.nvars = length(raw);
a.data = raw(:)';

% Set the var names.  These will be modified to make them valid, and the
% original strings saved in the VarDescription property.  Allow empty and
% duplicate names.
a = setvarnames(a,varNames(:)',[],true,true,true);


%-----------------------------------------------------------------------------
function raw = tdfread(file,delimiter,skip,treatAsEmpty)
%TDFREAD Read in text and numeric data from delimited file.

tab = sprintf('\t');
lf = sprintf('\n');
crlf = sprintf('\r\n');
cr = sprintf('\r');

% open file
fid = fopen(file,'rt'); % text mode: CRLF -> LF on windows (no-op on linux)
if fid == -1
    error(message('stats:dataset:dataset:OpenFailed', file));
end

% now read in the data
[bigM,count] = fread(fid,Inf);
fclose(fid);
if count == 0
    raw = cell(1,0);
    return
elseif bigM(count) ~= lf
   bigM = [bigM; lf];
end
bigM = char(bigM(:)');

% replace CRLF with LF (for reading DOS files on linux, where text mode is
% a no-op), then replace bare CR with LF (to handle macOS Excel files),
% then skip header lines, before removing empty lines below.
bigM = strrep(bigM,crlf,lf);
bigM = strrep(bigM,cr,lf);
if skip > 0
    i = find(bigM == lf,skip,'first');
    bigM(1:i(end)) = [];
end

% replace multiple embedded whitespace with a single whitespace, and multiple
% line breaks with one (removes empty lines in the middle and allows empty
% lines at the end).  remove insignificant whitespace before and after
% delimiters or line breaks.  remove leading empty line.
if delimiter == tab
%    matchexpr = {'([ \n])\1+' ' *(\n|\t) *' '^\n'};
   matchexpr = {'([\n])\1+' ' *(\n|\t) *' '^\n'};
elseif delimiter == ' '
%    matchexpr = {'([\t\n])\1+' '\t*(\n| )\t*' '^\n'};
   matchexpr = {'([\n])\1+' '\t*(\n| )\t*' '^\n'};
else
%    matchexpr = {'([\t\n])\1+' '\t*(\n| )\t*' '^\n'};
   matchexpr = {'([\n])\1+' ['[ \t]*(\n|\' delimiter ')[ \t]*'] '^\n'};
end
replexpr = {'$1' '$1' ''};
bigM = regexprep(bigM,matchexpr,replexpr);
delimitidx = find(bigM == delimiter);

% find out how many lines are there.
newlines = find(bigM == lf)';
nobs = length(newlines);

% find out how many data columns there are.
line1 = bigM(1:newlines(1)-1);
nvars = length(find(line1==delimiter)) + 1;

% check the size validation
if length(delimitidx) ~= nobs*(nvars-1)
   error(message('stats:dataset:dataset:BadFileFormat'));
end
if nvars > 1
   delimitidx = (reshape(delimitidx,nvars-1,nobs))';
end

startlines = [zeros(nobs>0,1); newlines(1:nobs-1)];
delimitidx = [startlines, delimitidx, newlines];
fieldlengths = diff(delimitidx,[],2) - 1; fieldlengths = fieldlengths(:);
if any(fieldlengths < 0)
   error(message('stats:dataset:dataset:BadFileFormat'));
end
maxlength = max(fieldlengths);
raw = cell(1,nvars);
for vars = 1:nvars
   xstr = repmat(' ',nobs,maxlength);
   x = NaN(nobs,1);
   xNumeric = true;
   for k = 1:nobs
       str = bigM(delimitidx(k,vars)+1:delimitidx(k,vars+1)-1);
       xstr(k,1:length(str)) = str;
       if xNumeric % numeric so far, anyways
           num = str2double(str);
           if isnan(num)
               % If the result was NaN, figure out why.  Note that because we
               % leave xstr alone, TreatAsEmpty only works on numeric columns.
               % That is what textscan does.
               if isempty(str) || any(strcmpi(str,'nan')) || any(strcmp(str,treatAsEmpty))
                   % Leave x(k) alone, it has a NaN already.  We won't decide
                   % if this column is numeric or not based on this value.
               else
                   % NaN must have come from a failed conversion, treat this
                   % column as strings.
                   xNumeric = false;
               end
           else
               % Otherwise accept the numeric value.  Numeric literals in
               % TreatAsEmpty, such as "-99", have already been disallowed, so
               % it's OK to accept any numeric literal here.  Numeric literals
               % cannot be used with TreatAsEmpty.  That is what textscan
               % does.
               x(k) = num;
           end
       end
   end
   if ~xNumeric
       x = cellstr(xstr);
   end
   raw{vars} = x;
end
