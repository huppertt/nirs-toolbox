function disp(t,bold,indent,fullChar)
%DISP Display a table.
%   DISP(T) prints the table T, including variable names and row names (if
%   present), without printing the table name.  In all other ways it's the
%   same as leaving the semicolon off an expression.
%
%   For numeric or categorical variables that are 2-dimensional and have 3 or
%   fewer columns, DISP prints the actual data using either short g, long g,
%   or bank format, depending on the current command line setting.  Otherwise,
%   DISP prints the size and type of each table element.
%
%   For character variables that are 2-dimensional and 10 or fewer characters
%   wide, DISP prints quoted strings.  Otherwise, DISP prints the size and
%   type of each table element.
%
%   For cell variables that are 2-dimensional and have 3 or fewer columns,
%   DISP prints the contents of each cell (or its size and type if too large).
%   Otherwise, DISP prints the size of each tabble element.
%
%   For other types of variables, DISP prints the size and type of each
%   table element.
%
%   See also TABLE, DISPLAY, FORMAT.

%   Copyright 2012-2013 The MathWorks, Inc.

if nargin < 2, bold = true; end
if nargin < 3, indent = 4; end
if nargin < 4, fullChar = false; end
between = 4;
within = 2;

% Follow the cmd window's format settings as possible
isLoose = strcmp(get(0,'FormatSpacing'),'loose');
if isLoose
    looseline = '\n';
else
    looseline = '';
end
[dblFmt,snglFmt] = getFloatFormats();

strongBegin = '';strongEnd = '';
if matlab.internal.display.isHot() && bold
    strongBegin = getString(message('MATLAB:table:localizedStrings:StrongBegin')); 
    strongEnd = getString(message('MATLAB:table:localizedStrings:StrongEnd')); 
end
varnameFmt = [strongBegin '%s' strongEnd];
    
if (t.nrows > 0) && (t.nvars > 0)
    indentSpaces = repmat(' ', t.nrows, indent);   % indent at left margin
    betweenSpaces = repmat(' ', t.nrows, between); % space betweeen variables
    withinSpaces = repmat(' ', t.nrows, within);   % space between columns within a variable
    if isempty(t.rownames)
        tblChars = indentSpaces;
    else
        rownameChars = char(t.rownames);
        rownameWidth = size(rownameChars,2);
        rownameChars = strcat(strongBegin,rownameChars,strongEnd);
        tblChars = [indentSpaces rownameChars betweenSpaces];
    end
    varnamePads = zeros(1,t.nvars);
    for ivar = 1:t.nvars
        name = t.varnames{ivar};
        var = t.data{ivar};
        
        if ischar(var)
            if ismatrix(var) && (fullChar || (size(var,2) <= 10))
                % Display individual strings for a char variable that is 2D and no
                % more than 10 chars.
                varChars = var;
            else
                % Otherwise, display a description of the chars.
                sz = size(var);
                szStr = ['[1' sprintf('x%d',sz(2:end)) ' char]'];
                varChars = repmat(szStr,sz(1),1);
            end
            
        else
            % Display the individual data if the var is 2D and no more than 3 columns.
            if ~isempty(var) && ismatrix(var) && (size(var,2) <= 3)
                if isnumeric(var)
                    if isa(var,'double')
                        varChars = num2str(var,dblFmt);
                    elseif isa(var,'single')
                        varChars = num2str(var,snglFmt);
                    else % integer types
                        varChars = num2str(var);
                    end
                elseif islogical(var)
                    % Display the logical values using meaningful names.
                    strs = ['false'; 'true '];
                    w1 = size(strs,2); w2 = within;
                    varChars = repmat(' ',size(var,1),(size(var,2)-1)*(w1+w2));
                    for j = 1:size(var,2)
                        varChars(:,(j-1)*(w1+w2)+(1:w1)) = strs(var(:,j)+1,:);
                    end
                elseif isa(var,'categorical') || isa(var,'datetime') || isa(var,'duration') || isa(var,'calendarDuration')
                    % Build the output one column at a time, since the char method reshapes
                    % to a single column.
                    varChars = char(zeros(t.nrows,0));
                    for j = 1:size(var,2)
                        if j > 1, varChars = [varChars withinSpaces]; end %#ok<AGROW>
                        varChars = [varChars char(var(:,j))]; %#ok<AGROW>
                    end
                elseif iscell(var)
                    % Let the built-in cell display method show the contents
                    % of each cell however it sees fit.  For example, it will
                    % display only a size/type if the contents are large.  It
                    % puts quotes around char contents, which char wouldn't.
                    varStr = evalc('disp(var)');
                    
                    % Work around a special case that the command line needs
                    % but we don't: curly braces around a scalar cell
                    % containing a 0x0
                    if isscalar(var) && max(size(var{1}))==0
                        varStr = removeBraces(varStr);
                    end
                    
                    % varStr is a single row with \n delimiting the chars for
                    % each row of var.  But \n can also be from displaying the
                    % contents of a cell.  There will be an extra trailing \n
                    % if isLoose; that can be left on.
                    loc = [0 find(varStr==10)];
                    [n,m] = size(var); % already checked is 2D
                    if length(loc) == n+1+isLoose % can use them as row delimiters
                        % The cell disp method puts leading whitespace
                        whiteSpace = find(varStr ~= ' ',1,'first') - 1;
                        % Split the \n-delimited string into a char matrix.
                        len = diff(loc) - whiteSpace;
                        varChars = repmat(' ',size(var,1),max(len)-1);
                        for i = 1:n
                            celChars = strtrim(varStr(loc(i)+1:loc(i+1)-1));
                            if ~isempty(celChars) % avoid 0x0 coming from strtrim
                                varChars(i,1:length(celChars)) = celChars;
                            end
                        end
                    else % the display for some cells had a \n in them
                        % Use the built-in to display each cell individually.
                        % This gives a slightly different output than the
                        % above, because cells are not all justified to the
                        % same length.
                        varChars = char(zeros(t.nrows,0));
                        offset = 0;
                        for j = 1:m
                            if j > 1
                                varChars = [varChars withinSpaces]; %#ok<AGROW>
                                offset = size(varChars,2);
                            end
                            for i = 1:n
                                % Display contents of each cell, remove {} around 0x0
                                var_ij = var(i,j);
                                celChars = evalc('disp(var_ij)');
                                if max(size(var_ij{1})) == 0
                                    celChars = removeBraces(celChars);
                                end
                                celChars = strtrim(celChars(1:end-1));
                                if ~isempty(celChars) % avoid 0x0 coming from strtrim
                                    varChars(i,offset+(1:length(celChars))) = celChars;
                                end
                            end
                        end
                    end
                elseif ~isempty(enumeration(var)) % isenumeration(var)
                    varChars = evalc('builtin(''disp'',var)');
                    if isLoose % remove trailing \n
                        varChars(end) = [];
                    end
                    numLines = size(var,1);
                    varChars = reshape(varChars,numel(varChars)/numLines,numLines)';
                    %Remove the name padding and trailing \n from enum DISP
                    varChars(:,[1:4 end]) = [];
                else
                    % Display a description of each table element.
                    varChars = getInfoDisplay(var);
                end

            % Either the variable is not 2D, or it's empty, or it's too wide
            % to show. Display a description of each table element.
            else
                varChars = getInfoDisplay(var);
            end
        end
        if size(varChars,2) < length(name)
            varChars(:,end+1:length(name)) = ' ';
        end
        varnamePads(ivar) = size(varChars,2)-length(name);
        
        if ivar == 1 % starting over at left margin
            tblChars = [tblChars varChars]; %#ok<AGROW>
        else
            tblChars = [tblChars betweenSpaces varChars]; %#ok<AGROW>
        end
    end
    dispVarNames();
    disp(tblChars);
else
    % do nothing
end
fprintf(looseline);

%-----------------------------------------------------------------------
    function dispVarNames()
        if isempty(t.rownames)
            fprintf('%s',repmat(' ',1,indent));
        else
            fprintf('%s',repmat(' ',1,indent+rownameWidth+between));
        end
        ii = 1;
        rightSpaces = repmat(' ',1,ceil(varnamePads(ii) / 2));
        leftSpaces = repmat(' ',1,varnamePads(ii)-length(rightSpaces));
        fprintf(varnameFmt,[leftSpaces t.varnames{ii} rightSpaces]);
        for ii = 2:t.nvars
            rightSpaces = repmat(' ',1,ceil(varnamePads(ii) / 2));
            leftSpaces = repmat(' ',1,varnamePads(ii)-length(rightSpaces)+between);
            fprintf(varnameFmt,[leftSpaces t.varnames{ii} rightSpaces]);
        end
        fprintf('\n');
        if isempty(t.rownames)
            fprintf('%s',repmat(' ',1,indent));
        else
            fprintf('%s',repmat(' ',1,indent+rownameWidth+between));
        end
        ii = 1;
        ul = repmat('_',1,length(t.varnames{ii})+varnamePads(ii));
        fprintf(varnameFmt,ul);
        for ii = 2:t.nvars
            spaces = repmat(' ',1,between);
            ul = repmat('_',1,length(t.varnames{ii})+varnamePads(ii));
            fprintf('%s',[spaces sprintf(varnameFmt,ul)]);
        end
        fprintf(['\n' looseline]);
    end

end % main function

%-----------------------------------------------------------------------
function [dblFmt,snglFmt] = getFloatFormats()
% Display for double/single will follow 'format long/short g/e' or 'format bank'
% from the command window. 'format long/short' (no 'g/e') is not supported
% because it often needs to print a leading scale factor.
switch get(0,'Format')
case {'short' 'shortG' 'shortEng'}
    dblFmt  = '%.5g    ';
    snglFmt = '%.5g    ';
case {'long' 'longG' 'longEng'}
    dblFmt  = '%.15g    ';
    snglFmt = '%.7g    ';
case 'shortE'
    dblFmt  = '%.4e    ';
    snglFmt = '%.4e    ';
case 'longE'
    dblFmt  = '%.14e    ';
    snglFmt = '%.6e    ';
case 'bank'
    dblFmt  = '%.2f    ';
    snglFmt = '%.2f    ';
otherwise % rat, hex, + fall back to shortg
    dblFmt  = '%.5g    ';
    snglFmt = '%.5g    ';
end
end


%-----------------------------------------------------------------------
function str = removeBraces(str)
str = regexprep(str,'\{(.*)\}','$1');
end


%-----------------------------------------------------------------------
function varChars = getInfoDisplay(var)
sz = size(var);
if ismatrix(var)
    szStr = ['[1' sprintf('x%d',sz(2:end))];
else
    szStr = ['[1' sprintf('x%d',sz(2)) sprintf('x%d',sz(3:end))];
end
varChars = repmat([szStr ' ' class(var) ']'],sz(1),1);
end
