function writeTextFile(t,file,args)
%WRITETEXTFILE Write a table to a text file.

%   Copyright 2012-2014 The MathWorks, Inc.

import matlab.internal.tableUtils.validateLogical

pnames = {'WriteVariableNames' 'WriteRowNames' 'Delimiter' 'QuoteStrings'};
dflts =  {               true          false       'comma'         false };
[writeVarNames,writeRowNames,delimiter,addQuotes] ...
                   = matlab.internal.table.parseArgs(pnames, dflts, args{:});

writeRowNames = validateLogical(writeRowNames,'WriteRowNames');
writeVarNames = validateLogical(writeVarNames,'WriteVariableNames');
addQuotes = validateLogical(addQuotes,'QuoteStrings');

% Only write row names if asked to, and if they exist.
writeRowNames = (writeRowNames && ~isempty(t.rownames));

tab = sprintf('\t');
lf = sprintf('\n');

% Set the delimiter, recognizing aliases
switch delimiter
case {'tab', '\t', tab}
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
  error(message('MATLAB:table:write:UnrecognizedDelimiter', delimiter(1)));
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
[fid,errmsg] = fopen(file,'wt'); % text mode: CRLF -> LF

if fid == -1
   error(message('MATLAB:table:write:FileOpenError', file, errmsg));
end

adata = t.data;
varType = zeros(1,t.nvars);
ncellCols = cell(1,t.nvars);
for j = 1:t.nvars
    varj = adata{j};
    varType(j) = typeCode(varj);
    if iscell(varj)
        ncellCols{j} = max(cellfun(@ncolsCell,varj),[],1);
    end
end

% Write variable names to the first line of the file as column headers
if writeVarNames
    str = '';

    % Start with a column header for row names
    if writeRowNames
        str = sprintf('%s%s',t.props.DimensionNames{1},delimiter);
    end

    for j = 1:t.nvars
        varnamej = colHeaders(adata{j},t.varnames(j),ncellCols{j});
        for jj = 1:length(varnamej)
            str = [str sprintf('%s%s',varnamej{jj},delimiter)]; %#ok<AGROW>
        end
    end

    % Write out the header line
    if ~isempty(str), str(end) = lf; end % replace trailing delimiter with '\n'
    fprintf(fid,'%s',str);
end

% Write each row of the table to the file
buf = ''; bufLen = 0;
for i = 1:t.nrows
    if writeRowNames
        str = sprintf('%s%s',t.rownames{i},delimiter);
        start = bufLen+1; bufLen = bufLen+length(str);
        buf(start:bufLen) = str;
    end
    for j = 1:t.nvars
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
if mod(t.nrows,10) ~= 0, fprintf(fid,'%s',buf(1:bufLen)); end

fclose(fid);


    function str = writeCell(x,ncellCols)
    % WRITECELL Write a cell-valued table element to the file
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
    % WRITEELEMENT Write a non-cell table element to the file
    
    if issparse(x), x = full(x); end % Only one row, no memory issues
    
    switch type
    case {1 2 3 4}
        if isreal(x)
            str = sprintf(realFormats{type},x); % may be multiple values
        else
            str = sprintf(complexFormats{type},[real(x); imag(x)]); % may be multiple values
        end
    case 5 % char
        if isrow(x) || isequal(x,'')
            if addQuotes
                str = ['"' strrep(x,'"','""') '"' delimiter];
            else
                str = [x delimiter];
            end
        else % multiple strings from a cell of a cellstr
            [n,m,d] = size(x);
            if d > 1 % permute/reshape N-D to matrix
                x = permute(x,[1 3:ndims(x) 2]);
                x = reshape(x,[n*d,m]);
            end
            % Turn them all into one long string.
            if addQuotes
                xc = cell(1,size(x,1));
                if isempty(xc)
                    str = delimiter; % output just a delimiter
                else
                    for s = 1:length(xc)
                        xc{s} = ['"' strrep(x(s,:),'"','""') '"' delimiter];
                    end
                    str = strjoin(xc,'');
                end
            else
                x(:,end+1) = delimiter;
                x = x'; str = x(:)';
            end
        end
    case {6 7}% categorical, datetime, duration, calendarDuration
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
function varnamej = colHeaders(varj,varnamej,ncellColsj)
%COLHEADERS Create multiple column headers from a table variable name

% Need multiple column headers if the variable has multiple columns.
if ischar(varj)
    [~,~,ncols] = size(varj);
else
    [~,ncols] = size(varj); % Treat N-D as 2-D.
end
if ncols > 1
    varnamej = strcat(varnamej,'_',num2str((1:ncols)'))';
end

% Need multiple column headers if the variable is a cell containing non-scalars.
if iscell(varj) && any(ncellColsj(:) > 1)
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
elseif isa(x,'datetime') || isa(x,'duration') || isa(x,'calendarDuration')
    type = 7;
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
