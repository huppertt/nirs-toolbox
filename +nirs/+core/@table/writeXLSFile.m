function writeXLSFile(t,xptfile,args)
%WRITEXLSFILE Write a table to an Excel spreadsheet file.

%   Copyright 2012-2014 The MathWorks, Inc.

import matlab.internal.tableUtils.validateLogical
import matlab.internal.tableUtils.matricize

pnames = {'WriteVariableNames' 'WriteRowNames' 'Sheet' 'Range'};
dflts =  {                true          false        1     'A1'};
[writeVarNames,writeRowNames,sheet,range] ...
                   = matlab.internal.table.parseArgs(pnames, dflts, args{:});

writeRowNames = validateLogical(writeRowNames,'WriteRowNames');
writeVarNames = validateLogical(writeVarNames,'WriteVariableNames');

% Only write row names if asked to, and if they exist.
writeRowNames = (writeRowNames && ~isempty(t.rownames));

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
   error(message('MATLAB:table:write:InvalidRange'));
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

% Write row names.
vars = cell(t.nrows+writeVarNames,0);
ulcol = col;
if writeRowNames
    rownames = t.rownames;
    if writeVarNames
        rownames = [t.props.DimensionNames{1}; rownames];
    end
    vars = [vars rownames];
    col = col + 1;
end

for j = 1:t.nvars
    if col > maxCol, break; end
    
    varj = t.data{j};
    varnamej = t.varnames{j};

    if iscell(varj)
        % xlswrite cannot write out non-scalar-valued cells -- convert cell
        % variables to a cell of the appropriate width containing only
        % scalars.
        [~,ncols] = size(varj); % Treat N-D as 2-D
        ncellColsj = max(cellfun(@ncolsCell,matricize(varj)),[],1);
        newNumCols = sum(ncellColsj);
        newVarj = cell(t.nrows,newNumCols);
        
        % Expand out each column of varj into as many columns as needed to
        % have only scalar-valued cells, possibly padded with empty cells.
        cnt = 0;
        for jj = 1:ncols
            varjj = varj(:,jj);
            num = ncellColsj(jj);
            newVarjj = cell(t.nrows,num);
            for i = 1:t.nrows
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
                elseif isa(varjj_i,'categorical') || isa(varjj_i,'datetime') ||  isa(varjj_i,'duration') ||  isa(varjj_i,'calendarDuration')
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
            varj = num2cell(matricize(varj));
        elseif isa(varj,'categorical') || isa(varj,'datetime') ||  isa(varj,'duration') ||  isa(varj,'calendarDuration')
            varj = cellstr(matricize(varj));
        else % write out empty cells
            varj = cell(t.nrows,1);
        end
    end

    [~,ncols] = size(varj); % Treat N-D as 2-D
    if writeVarNames
        if ncols > 1
            varj = [strcat({varnamej},'_',num2str((1:ncols)'))'; varj]; %#ok<AGROW>
        else
            varj = [{varnamej}; varj]; %#ok<AGROW>
        end
    end
    vars = [vars varj]; %#ok<AGROW>
    col = col + ncols;

    if mod(j,10) == 0
        writeXLSVars();
        vars = cell(t.nrows+writeVarNames,0);
        ulcol = col;
    end
end
if mod(j,10) ~= 0, writeXLSVars(); end

    function writeXLSVars
        ul = rowcol2spec(row,ulcol);
        [success,msg] = xlswrite(xptfile,vars,sheet,[ul lr]);
        if ~success
            error(message('MATLAB:table:write:XlsWriteFailed', varnamej, xptfile, msg.message));
        end
    end

end % writeXLSFile function


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
    if ~ismatrix(c)
        c = permute(c,[2 1 3:ndims(c)]);
        c = reshape(c,size(c,1),[])';
    end
    cs = reshape(cellstr(c),szOut);
end
end


%-----------------------------------------------------------------------
function spec = rowcol2spec(row,col)
col = min(col,16384); % 'xfd'
mod26 = mod(col-1,26);
if col <= 26
    spec = sprintf('%s%d',char('A'+mod26),row);
else
    mult26 = mod(floor((col-26-1)/26),26);
    if col <= 702
        spec = sprintf('%s%s%d',char('A'+mult26),char('A'+mod26),row);
    else
        mult676 = floor((col-702-1)/676);
        spec = sprintf('%s%s%s%d',char('A'+mult676),char('A'+mult26),char('A'+mod26),row);
    end
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
