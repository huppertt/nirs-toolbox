function export(a,varargin)
%EXPORT Write a dataset array to a file.
%   EXPORT(DS,'File',FILENAME) writes the dataset array DS to a tab-delimited
%   text file, including variable names and observation names (if present).
%   If the observation names exist, the name in the first column of the first
%   line of the file is the first dimension name for the dataset (by default,
%   'Observations').  EXPORT overwrites any existing file named FILENAME.
%
%   EXPORT(DS) writes to a text file whose default name is the name of the
%   dataset array DS, appended by '.txt'.  If EXPORT cannot construct the file
%   name from the dataset array input, it writes to the file 'dataset.txt'.
%   EXPORT overwrites any existing file.
%
%   EXPORT(DS,'File',FILENAME,'Delimiter',DELIM) writes the dataset array DS
%   to a text file, using the delimiter DELIM.  DELIM can be any of ' ', '\t',
%   ',', ';', '|' or their corresponding string names 'space', 'tab', 'comma',
%   'semi', or 'bar'.
%
%   EXPORT(DS,'XLSFile',FILENAME) writes the dataset array DS to an Excel
%   spreadsheet file, including variable names and observation names (if
%   present).  You may also specify the 'Sheet' and 'Range' parameter
%   name/value pairs, with parameter values as accepted by the XLSREAD
%   function.
%
%   EXPORT(DS,'XPTFile',FILENAME) writes the dataset array DS to a SAS
%   XPORT format file. When writing to an XPORT format file variables must
%   be scalar valued. EXPORT saves observation names to a variable called
%   OBSNAMES unless the WriteObsNames parameter described below is set to
%   false. The XPORT format restricts the length of variable names to eight
%   characters; longer variable names will be truncated.
%
%   EXPORT(DS, ..., 'WriteVarNames',FALSE) does not write the variable names
%   to the text file.  EXPORT(..., 'WriteVarNames',TRUE) writes the names as
%   column headings in the first line of the file, and is the default.
%
%   EXPORT(DS ..., 'WriteObsNames',FALSE) does not write the observation names
%   to the text file.  EXPORT(..., 'WriteObsNames',TRUE) writes the names as
%   the first column of the file, and is the default.
%
%   In some cases, EXPORT creates a text file that does not represent DS
%   exactly, as described below.  If you use DATASET('File',FILENAME) to read
%   that file back in and create a new dataset array, the result may not have
%   exactly the same format or contents as the original dataset array.
%
%   *  EXPORT writes out numeric variables using long g format, and
%      categorical or character variables as unquoted strings.
%   *  For non-character variables that have more than one column, EXPORT
%      writes out multiple delimiter-separated fields on each line, and
%      constructs suitable column headings for the first line of the file.
%   *  EXPORT writes out variables that have more than two dimensions as two
%      dimensional variables, with trailing dimensions collapsed.
%   *  For cell-valued variables, EXPORT writes out the contents of each cell
%      as a single row, in multiple delimiter-separated fields, when the
%      contents are numeric, logical, character, or categorical, and writes
%      out a single empty field otherwise.
%
%   Save DS as a mat file if you need to import it again as a dataset array.
%
%   Examples:
%      load hospital
%
%      % Write to a comma-delimited file 'hospital.txt'
%      export(hospital,'delimiter',',')
%
%      % Write to a comma-delimited file 'dataset.txt'
%      export(hospital(1:10,1:3),'delimiter',',')
%
%      % Write to a comma-delimited file 'Under40.dat'
%      export(hospital(hospital.Age<40,:),'File','Under40.dat','delimiter',',')
%      
%   See also DATASET.

%   Copyright 2008-2012 The MathWorks, Inc.

%   SAS is a registered trademarks of SAS Institute Inc.

pnames = {'file' 'xlsfile' 'xptfile' 'writevarnames' 'writeobsnames' };
dflts =  {   []        []      []      true            true  };
[fileArg,xlsfileArg,xptfileArg,writevarnames,writeobsnames,supplied,otherArgs] = ...
    dataset.parseArgs(pnames, dflts, varargin{:});

if ~isscalar(writevarnames) || ...
                ~(islogical(writevarnames) || isnumeric(writevarnames))
    error(message('stats:dataset:export:InvalidWriteVarNames'));
end
if ~isscalar(writeobsnames) || ...
                ~(islogical(writeobsnames) || isnumeric(writeobsnames))
    error(message('stats:dataset:export:InvalidWriteObsNames'));
end

% Only write obsnames if asked to, and if they exist.
writeobsnames = (writeobsnames && ~isempty(a.obsnames));

if supplied.file || (~supplied.xlsfile && ~supplied.xptfile)
    % Create a default name if needed
    if ~supplied.file
        dsname = inputname(1);
        if isempty(dsname)
            dsname = 'dataset';
        end
        fileArg = [dsname '.txt'];
    end
    writeFile(a,fileArg,writevarnames,writeobsnames,otherArgs{:});
elseif supplied.xlsfile
    try
        Excel = actxserver('Excel.Application');
    catch me %#ok<NASGU>
        error(message('stats:dataset:export:NoCOMServer'));
    end
    Excel.Quit;
    
    writeXLSFile(a,xlsfileArg,writevarnames,writeobsnames,otherArgs{:});
else % Write to xpt file
    xptwrite(a,xptfileArg,'WriteObsNames',writeobsnames,...
             'DataSetName',inputname(1),otherArgs{:});
end

end


%-----------------------------------------------------------------------
function writeFile(a,filename,writevarnames,writeobsnames,varargin)
%WRITEFILE Write a dataset array to a text file

pnames = {'delimiter'};
dflts =  {      'tab'};
delimiter = dataset.parseArgs(pnames, dflts, varargin{:});

tab = sprintf('\t');
lf = sprintf('\n');

% Set the delimiter
switch delimiter
case {'tab', '\t'}
  delimiter = tab;
case {'space',' '}
  delimiter = ' ';
case {'comma', ','}
  delimiter = ',';
case {'semi', ';'}
  delimiter = ';';
case {'bar', '|'}
  delimiter = '|';
otherwise
  error(message('stats:dataset:export:UnrecognizedDelimiter', delimiter( 1 )));
end

realDoubleFmt = ['%.15g' delimiter];
complexDoubleFmt = ['%.15g+%.15gi' delimiter];
realSingleFmt = ['%.7g' delimiter];
complexSingleFmt = ['%.7g+%.7gi' delimiter];
realIntegerFmt = ['%d' delimiter];
complexIntegerFmt = ['%d+%di' delimiter];
realFormats = {realDoubleFmt realSingleFmt realIntegerFmt realIntegerFmt};
complexFormats = {complexDoubleFmt complexSingleFmt complexIntegerFmt complexIntegerFmt};

% Open the file for writing.
fid = fopen(filename,'wt'); % text mode: CRLF -> LF

if fid == -1
   error(message('stats:dataset:export:ExportOpenFailed', filename));
end

adata = a.data;
varType = zeros(1,a.nvars);
ncellCols = cell(1,a.nvars);
for j = 1:a.nvars
    varj = adata{j};
    varType(j) = typeCode(varj);
    if iscell(varj)
        ncellCols{j} = max(cellfun(@ncolsCell,varj),[],1);
    end
end

% Write variable names to the first line of the file as column headers
if writevarnames
    str = '';

    % Start with a column header for observation names
    if writeobsnames
        str = sprintf('%s%s',a.props.DimNames{1},delimiter);
    end

    for j = 1:a.nvars
        varnamej = colHeaders(adata{j},a.varnames(j),ncellCols{j});
        for jj = 1:length(varnamej)
            str = [str sprintf('%s%s',varnamej{jj},delimiter)];
        end
    end

    % Write out the header line
    if ~isempty(str), str(end) = lf; end % replace trailing delimiter with '\n'
    fprintf(fid,'%s',str);
end

% Write each row of the dataset array to the file
buf = ''; bufLen = 0;
for i = 1:a.nobs
    if writeobsnames
        str = sprintf('%s%s',a.obsnames{i},delimiter);
        start = bufLen+1; bufLen = bufLen+length(str);
        buf(start:bufLen) = str;
    end
    for j = 1:a.nvars
        varj = adata{j};
        type = varType(j);
        if type == 8 % cell variable
            str = writeCell(varj(i,:),ncellCols{j});
        elseif type == 5 % char variable
            str = writeElement(varj(i,:,:),type); % keep each row separate
        else % type < 8, non-cell variable
            str = writeElement(varj(i,:),type); % N-D to 2-D
        end
        start = bufLen+1; bufLen = bufLen+length(str);
        buf(start:bufLen) = str;
    end
    if ~isempty(buf), buf(bufLen) = lf; end % replace trailing delimiter with '\n'
    if mod(i,10) == 0
        fprintf(fid,'%s',buf(1:bufLen));
        bufLen = 0;
    end
end
if mod(a.nobs,10) ~= 0, fprintf(fid,'%s',buf(1:bufLen)); end

fclose(fid);


    function str = writeCell(x,ncellCols)
    % WRITECELL Write a cell-valued dataset array element to the file
        [~,ncols] = size(x); % Treat N-D as 2-D
        str = '';
        for k = 1:ncols
            cellk = x{k};
            cellType = typeCode(x{k});
            
            % Pad with empty fields if necessary to match other rows.
            pad = delimiter(1,ones(1,ncellCols(k)-ncolsCell(cellk)));
            if cellType == 5
                % Write the cell as a set of char strings.  
                str = sprintf('%s%s%s',str,writeElement(cellk,cellType),pad);
            else
                % Write the cell as a row of values.  writeElement will write an
                % empty field if the contents are not atomic.
                str = sprintf('%s%s%s',str,writeElement(cellk(:)',cellType),pad);
            end
        end
    end

    function str = writeElement(x,type)
    % WRITEELEMENT Write a non-cell dataset array element to the file
    
    if issparse(x), x = full(x); end % Only one row, no memory issues
    
    switch type
    case {1 2 3 4}
        if isreal(x)
            str = sprintf(realFormats{type},x); % may be multiple values
        else
            str = sprintf(complexFormats{type},[real(x); imag(x)]); % may be multiple values
        end
    case 5 % char
        [n,~,d] = size(x);
        if n*d == 1 % a single string
            str = [x,delimiter];
        elseif n*d == 0 % '' or zero strings
            str = delimiter;
        else % n*d > 1: multiple strings, including N-D
            % Write out each row as a separate field in the string, including
            % rows in higher dims.
            x = permute(x,[2 1 3:ndims(x)]);
            x = x(:,:); x(end+1,:) = delimiter;
            str = x(:)';
        end
    case 6 % categorical
        if size(x,2) == 1
            str = [char(x) delimiter];
        else
            % Can't use strcat, that drops trailing whitespace, including tabs
            str = cell2mat(cellfun(@(str)[str delimiter],cellstr(x),'UniformOutput',false));  % may be multiple values
        end
    otherwise % not a known type, or a cell
        % Don't attempt to write out
        str = sprintf('%s',delimiter);
    end
    end

end % writeFile function


%-----------------------------------------------------------------------
function writeXLSFile(a,filename,writevarnames,writeobsnames,varargin)
%WRITEXLSFILE Write a dataset array to an Excel spreadsheet file

pnames = {'sheet' 'range'};
dflts =  {     1     'A1'};
[sheet,range] = dataset.parseArgs(pnames, dflts, varargin{:});

lr = '';
validRange = false;
if ischar(range) && isvector(range) && size(range,1)==1
    rangeSplit = regexpi(range,':','split');
    if length(rangeSplit)==1 || length(rangeSplit)==2
        ul = rangeSplit{1};
        [row,col] = spec2rowcol(ul);
        validRange = ~(isnan(row) || isnan(col));
        if validRange && length(rangeSplit) > 1
            lr = rangeSplit{2};
            [row2,col2] = spec2rowcol(lr);
            validRange = ~(isnan(row2) || isnan(col2));
        end
    end
end
if ~validRange
   error(message('stats:dataset:export:InvalidRange'));
end
if isempty(lr)
    maxCol = Inf;
else
    % convert ur:ll to ul:lr
    tmp = min(col,col2); col2 = max(col,col2); col = tmp;
    tmp = min(row,row2); row2 = max(row,row2); row = tmp;
    lr = [':' rowcol2spec(row2,col2)];
    maxCol = col2;
end

% Write observation names.
vars = cell(a.nobs+writevarnames,0);
ulcol = col;
if writeobsnames
    obsnames = a.obsnames;
    if writevarnames
        obsnames = [a.props.DimNames{1}; obsnames];
    end
    vars = [vars obsnames];
    col = col + 1;
end

for j = 1:a.nvars
    if col > maxCol, break; end
    
    varj = a.data{j};
    varnamej = a.varnames(j);

    if iscell(varj)
        % xlswrite cannot write out non-scalar-valued cells -- convert cell
        % variables to a cell of the appropriate width containing only
        % scalars.
        [~,ncols] = size(varj); % Treat N-D as 2-D
        ncellColsj = max(cellfun(@ncolsCell,varj),[],1);
        newNumCols = sum(ncellColsj);
        newVarj = cell(a.nobs,newNumCols);
        
        % Expand out each column of varj into as many columns as needed to
        % have only scalar-valued cells, possibly padded with empty cells.
        cnt = 0;
        for jj = 1:ncols
            varjj = varj(:,jj);
            num = ncellColsj(jj);
            newVarjj = cell(a.nobs,num);
            for i = 1:a.nobs
                % Expand each cell with non-scalar contents into a row of cells containing scalars
                varjj_i = varjj{i};
                if ischar(varjj_i)
                    % Put each string into its own cell.  If there are no
                    % strings (zero rows or zero pages in the original char
                    % array), the output will be a single empty cell.
                    vals = char2cell(varjj_i); % creates a 2-D cellstr
                    if isempty(vals), vals = {''}; end
                elseif isnumeric(varjj_i) || islogical(varjj_i)
                    vals = num2cell(varjj_i);
                elseif isa(varjj_i,'categorical')
                    vals = cellstr(varjj_i);
                else
                    vals = cell(0,0); % write out only an empty cell
                end
                newVarjj(i,1:numel(vals)) = vals(:)';
            end
            newVarj(:,cnt+(1:num)) = newVarjj;
            cnt = cnt + num;
        end
        
        varj = newVarj;
    else
        % xlswrite will convert any input to cell array anyway, may as well do
        % it here in all cases to get correct behavior for character and for
        % cases xlswrite won't handle.
        if ischar(varj)
            varj = char2cell(varj);
        elseif isnumeric(varj) || islogical(varj)
            varj = num2cell(varj(:,:));
        elseif isa(varj,'categorical')
            varj = cellstr(varj(:,:));
        else % write out empty cells
            varj = cell(a.nobs,1);
        end
    end

    [~,ncols] = size(varj); % Treat N-D as 2-D
    if writevarnames
        if ncols > 1
            varnamej = strcat(varnamej,'_',num2str((1:ncols)'))';
        end
        varj = [varnamej; varj];
    end
    vars = [vars varj];
    col = col + ncols;

    if mod(j,10) == 0
        writeXLSVars();
        vars = cell(a.nobs+writevarnames,0);
        ulcol = col;
    end
end
if mod(j,10) ~= 0, writeXLSVars(); end

    function writeXLSVars
        ul = rowcol2spec(row,ulcol);
        [success,msg] = xlswrite(filename,vars,sheet,[ul lr]);
        if ~success
            error(message('stats:dataset:export:XlswriteFailed', a.varnames{ j }, filename, msg.message));
        end
    end

end % writeXLSFile function


%-----------------------------------------------------------------------
function varnamej = colHeaders(varj,varnamej,ncellColsj)
%COLHEADERS Create multiple column headers from a dataset variable name

% Need multiple column headers if the variable has multiple columns.
if ischar(varj)
    ncols = 1;
else
    [~,ncols] = size(varj); % Treat N-D as 2-D.
end
if ncols > 1
    varnamej = strcat(varnamej,'_',num2str((1:ncols)'))';
end

% Need multiple column headers if the variable is a cell containing non-scalars.
if iscell(varj) && any(ncellColsj > 1)
    vnj = cell(1,sum(ncellColsj));
    cnt = 0;
    for jj = 1:ncols
        num = ncellColsj(jj);
        vnj(cnt+(1:num)) = strcat(varnamej(jj),'_',num2str((1:num)'))';
        cnt = cnt + num;
    end
    varnamej = vnj;
end

end % colHeaders function


%-----------------------------------------------------------------------
function type = typeCode(x)
% TYPECODE Return the type of a variable, coded as an integer.
if isa(x,'double')
    type = 1;
elseif isa(x,'single')
    type = 2;
elseif isinteger(x)
    type = 3;
elseif islogical(x)
    type = 4;
elseif ischar(x)
    type = 5;
elseif isa(x,'categorical')
    type = 6;
% elseif isa(x,'timeseries') % special handling of old-style object has been removed
%     type = 7;
elseif iscell(x)
    type = 8;
else
    type = 0; % other, in this case not a standard type
end
end


%-----------------------------------------------------------------------
function m = ncolsCell(c)
% How many columns will be needed to write out the contents of a cell?
if ischar(c)
    % Treat each row as a separate string, including rows in higher dims.
    [n,~,d] = size(c);
    % Each string gets one "column".  Zero rows (no strings) gets a single
    % column to contain the empty string, even for N-D,.  In particular,
    % '' gets one column.
    m = max(n*d,1);
elseif isnumeric(c) || islogical(c) || isa(c,'categorical')
    m = max(numel(c),1); % always write out at least one empty field
else
    m = 1; % other types are written as an empty field
end
end


%-----------------------------------------------------------------------
function cs = char2cell(c)
% Convert a char array to a cell array of strings, each cell containing a
% single string.  Treat each row as a separate string, including rows in
% higher dims.

% Create a cellstr array the same size as the original char array (ignoring
% columns), except with trailing dimensions collapsed down to 2-D.
[n,~,d] = size(c); szOut = [n,d];

if isempty(c)
    % cellstr would converts any empty char to {''}.  Instead, preserve the
    % desired size.
    cs = repmat({''},szOut);
else
    % cellstr does not accept N-D char arrays, put pages as more rows.
    if ndims(c) > 2
        c = permute(c,[2 1 3:ndims(c)]);
        c = c(:,:)';
    end
    cs = reshape(cellstr(c),szOut);
end
end


%-----------------------------------------------------------------------
function spec = rowcol2spec(row,col)
mult676 = floor((col-1)/676);
mult26 = floor((col-mult676*676-1)/26);
mod26 = mod(col-1,26);
if col <= 26
    spec = sprintf('%s%d',char('A'+mod26),row);
elseif col <= 702
    spec = sprintf('%s%s%d',char('A'+mult26-1),char('A'+mod26),row);
else
    spec = sprintf('%s%s%s%d',char('A'+mult676-1),char('A'+mult26-1),char('A'+mod26),row);
end
end


%-----------------------------------------------------------------------
function [row,col] = spec2rowcol(corner)
row = NaN; col = NaN;
tok = regexpi(corner,'^([a-z]{1,3})([0-9]+)$','tokens');
if ~isempty(tok)
    tok = tok{1};
    if length(tok) == 2
        c = tok{1}; r = tok{2};
        mults = [676; 26; 1]; % 'a' => 1, 'xfd' => 16384
        col = (lower(c) - 'a' + 1) * mults((end-length(c)+1):end);
        row = str2double(r);
    end
    % Let xlswrite decide if the row or column is too large
end
end
