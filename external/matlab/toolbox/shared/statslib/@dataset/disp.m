function disp(a)
%DISP Display a dataset array.
%   DISP(DS) prints the dataset array DS, including variable names and
%   observation names (if present), without printing the dataset name.  In all
%   other ways it's the same as leaving the semicolon off an expression.
%
%   For numeric or categorical variables that are 2-dimensional and have 3 or
%   fewer columns, DISP prints the actual data using either short g, long g,
%   or bank format, depending on the current command line setting.  Otherwise,
%   DISP prints the size and type of each dataset element.
%
%   For character variables that are 2-dimensional and 10 or fewer characters
%   wide, DISP prints quoted strings.  Otherwise, DISP prints the size and
%   type of each dataset element.
%
%   For cell variables that are 2-dimensional and have 3 or fewer columns,
%   DISP prints the contents of each cell (or its size and type if too large).
%   Otherwise, DISP prints the size of each dataset element.
%
%   For other types of variables, DISP prints the size and type of each
%   dataset element.
%
%   See also DATASET, @DATASET/DISPLAY, FORMAT.

%   Copyright 2006-2013 The MathWorks, Inc. 


isLoose = strcmp(get(0,'FormatSpacing'),'loose');

% Display for double/single will approximate 'format long/short g' or 'format
% bank'. 'format long/short e' might also be reasonable.  'format long/short'
% (no 'g') is not supported because that often needs to print a leading scale
% factor.
isBank = strcmp(get(0,'Format'),'bank');
if isBank
    dblFmt = '%.2f    ';
    snglFmt = '%.2f    ';
else
    isLong = ~isempty(strfind(get(0,'Format'),'long'));
    dblFmt = 5 + 10*isLong; % 5 or 15
    snglFmt = 5 + 2*isLong; % 5 or 7
end

maxWidth = matlab.desktop.commandwindow.size; maxWidth = maxWidth(1);

if (a.nobs > 0) && (a.nvars > 0)
    dsPad = repmat(' ', a.nobs+1, 4);
    varPad = repmat(' ', a.nobs, 2);
    if isempty(a.obsnames)
        dsChars = char(zeros(a.nobs+1, 0));
    else
        dsChars = [dsPad char(' ',char(a.obsnames))];
    end
    for ivar = 1:a.nvars
        name = a.varnames{ivar};
        var = a.data{ivar};
        
        if ischar(var)
            if (ndims(var) == 2) && (size(var,2) <= 10)
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
            if ~isempty(var) && (ndims(var) == 2) && (size(var,2) <= 3)
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
                    w1 = size(strs,2); w2 = size(varPad,2);
                    varChars = repmat(' ',size(var,1),(size(var,2)-1)*(w1+w2));
                    for j = 1:size(var,2)
                        varChars(:,(j-1)*(w1+w2)+(1:w1)) = strs(var(:,j)+1,:);
                    end
                elseif isa(var,'categorical')
                    % Build the output one column at a time, since the char method reshapes
                    % to a single column.
                    varChars = char(zeros(a.nobs,0));
                    for j = 1:size(var,2)
                        if j > 1, varChars = [varChars varPad]; end
                        varChars = [varChars char(var(:,j))];
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
                        % Split the \n-delimited string into a char matrix.
                        len = diff(loc);
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
                        varChars = char(zeros(a.nobs,0));
                        offset = 0;
                        for j = 1:m
                            if j > 1
                                varChars = [varChars varPad];
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
                else
                    % Display a description of each dataset element.
                    sz = size(var);
                    szStr = ['[1' sprintf('x%d',sz(2:end))];
                    varChars = repmat([szStr ' ' class(var) ']'],sz(1),1);
                end

            % Either the variable is not 2D, or it's empty, or it's too wide
            % to show. Display a description of each dataset element.
            else
                sz = size(var);
                if ndims(var) == 2
                    szStr = ['[1' sprintf('x%d',sz(2:end))];
                else
                    szStr = ['[' sprintf('%d',sz(2)) sprintf('x%d',sz(3:end))];
                end
                varChars = repmat([szStr ' ' class(var) ']'],sz(1),1);
            end
        end
        varChars = char(name,varChars);
        
        % If this new variable will extend the display past the right margin
        % with, display the output built up so far, and then restart for
        % display starting at the left margin.  Don't do that if this is the
        % first variable, otherwise we'd display only the observation names.
        if ivar > 1 && size(dsChars,2) + size(dsPad,2) + size(varChars,2) > maxWidth
                disp(dsChars);
                fprintf('\n');
                if (isLoose), fprintf('\n'); end
                if isempty(a.obsnames)
                    dsChars = char(zeros(a.nobs+1, 0));
                else
                    dsChars = [dsPad char(' ',a.obsnames{:})];
                end
        end
        dsChars = [dsChars dsPad varChars];
    end
    disp(dsChars);
else
    % do nothing
end
if (isLoose), fprintf('\n'); end


%-----------------------------------------------------------------------
function str = removeBraces(str)
str = regexprep(str,'\{(.*)\}','$1');
