function [str,msg] = variableEditorSetDataCode(a,varname,row,col,rhs)
%   This function is undocumented and will change in a future release

% Generate MATLAB command to edit the content of a cell to the specified
% rhs.

%   Copyright 2011-2014 The MathWorks, Inc.

msg = '';
str = '';

[colNames,varIndices,colClasses] = variableEditorColumnNames(a);

if length(row) > 1
    row = row(1);
end

singleCell = (length(col) == 1);
if length(col) > 1
    col = col(1);
end

if col>varIndices(end)
    msg = getString(message('MATLAB:table:VarEditorIndexOverflow'));
    col = varIndices(end);
end

% Find the left hand column of the table variable that contains col
j = find(varIndices>=col+1,1,'first')-1;

if ~isempty(j)
    colname = colNames{j};
    vardata = a.data{j};
    isGroupColumn = varIndices(j+1)>varIndices(j)+1;
    isCatColumn = isa(vardata,'categorical');
    emptyRHS = isempty(rhs);
    if isCatColumn
        % If the new categorical is just an empty string (''), set it to
        % <undefined>
        if strcmp(rhs, '')
            rhs = '<undefined>';
        end

        % Wrap the RHS string in quotes so that it appears as a string literal
        % in the generated code. Double any quotes in the RHS string so they act
        % as single quotes when rhs is itself put between quotes.
        rhs = ['''' strrep(rhs,'''','''''') ''''];
    end
    
    % If rhs is empty the evalin will error out with a 
    % "One or more output arguments not assigned" error which is not in
    % context for the user. Prevent this by replacing empty by the
    % appropriate empty designation.
    if isempty(rhs)
        if strcmp(colClasses{j},'cell')
            rhs = '{}';
        elseif strcmp(colClasses{j}, 'double') || ...
                strncmp(colClasses{j}, 'uint', 4) || ...
                strncmp(colClasses{j}, 'int', 3)
            rhs = '0';
        else
            rhs = '[]';
        end
    end
    rhsValue = evalin('caller',rhs);
    
    % Generate code only when rhs is a new value. This avoids cluttering the
    % command window when the table element would not change. Upstream code
    % catches a rhs that is literally the same, here we catch cases where it is
    % an expression that evaluates to the same value.
        
    % If this is a grouped column, an additional numeric sub-index is needed
    % to define which column within the table is being edited. A char matrix
    % is treated as one column, not as grouped columns.
    if isGroupColumn
        relcol = col-varIndices(j)+1;
        if iscell(vardata) % e.g. x.BloodPressure{10,2} = rhs;
            % Generate code only if the new value is different
            subsrefStruct = struct('type',{'.','{}'},'subs',{colname,{row,relcol}});
            if emptyRHS && ~singleCell
                % can't assign mutiple values to some table variable types
                % with a single value, so build a command with multiple
                % empties:  t.column(1,:) = {{},{}};
                str = [varname '.' colname '(' num2str(row) ',:) = {'];
                numEmpties = varIndices(j+1) - varIndices(j);
                for i=1:numEmpties
                    str = [str '{},']; %#ok<AGROW>
                end
                str = [str(1:length(str)-1) '};'];
            elseif row>size(a,1) || ~isequal(subsref(a,subsrefStruct),rhsValue)
                str = [varname '.' colname '{' num2str(row) ',' num2str(relcol) '} = ' rhs ';'];
            end
        else % e.g. x.BloodPressure(10,2) = rhs;
            if isa(vardata,'categorical') 
                % Already enforced that rhsValue is a string for a categorical var
            elseif ~isscalar(rhsValue)
                if isempty(rhsValue)
                    error(message('MATLAB:datatypes:InvalidValue', colname));
                else
                    error(message('MATLAB:subsassignnumelmismatch'));
                end
            end
            % Generate code only if the new value is different
            subsrefStruct = struct('type',{'.','()'},'subs',{colname,{row,relcol}});
            if emptyRHS && ~singleCell
                str = [varname '.' colname '(' num2str(row) ',:) = ' rhs ';'];
            elseif row>size(a,1) || ~isequal(subsref(a,subsrefStruct),rhsValue)
                str = [varname '.' colname '(' num2str(row) ',' num2str(relcol) ') = ' rhs ';'];
            end
        end
    else
        if iscell(vardata) % e.g. x.Weight{10} = rhs;
            % Generate code only if the new value is different
            subsrefStruct = struct('type',{'.','{}'},'subs',{colname,{row}});
            if row>size(a,1) || ~isequal(subsref(a,subsrefStruct),rhsValue)
                if size(a,1)==1 && row>1 % force the var to grow vertically
                    str = [varname '.' colname '{' num2str(row) ',1} = ' rhs ';'];
                else
                    str = [varname '.' colname '{' num2str(row) '} = ' rhs ';'];
                end
            end
        elseif strcmp('char',colClasses{j})
            % Need to get all the content of the value at the specified row
            % of the character array in order to compare against the new
            % value.
            if isempty(rhsValue)
                error(message('MATLAB:datatypes:InvalidValue', colname));
            else
                [~, columns] = size(subsref(a, struct('type',{'.'},'subs',{colname})));
                if columns == 1
                    subsrefStruct = struct('type',{'.','()'},'subs',{colname,{row, 1}});
                else
                    subsrefStruct = struct('type',{'.','()'},'subs',{colname,{row, [1, columns]}});
                end
                if row>size(a,1) || ~isequal(subsref(a,subsrefStruct),rhsValue)
                    str =  [varname '.' colname '(' num2str(row) ',1:' num2str(length(rhsValue)) ') = ' rhs ';'];
                end
            end
        elseif strcmp('datetime', colClasses{j})
            rhs = strtrim(rhs);
            if ismember(lower(rhs),{'''now''' '''yesterday''' '''today''' '''tomorrow'''}) || ...
                ismember(lower(rhs),{'now' 'yesterday' 'today' 'tomorrow'})
                % Allow keyword entry like the datetime constructor does
                str = [varname '.' colname '(' num2str(row) ') = datetime(' rhs ');'];
            elseif isempty(rhsValue)
                error(message('MATLAB:datatypes:InvalidValue', colname));
            else
                subsrefStruct = struct('type',{'.','()'},'subs',{colname,{row}});

                % Verify this is a valid datetime, with the format of the
                % current datetime array.  If it is, assign using the
                % string - the datetime constructor will use the current
                % format by default.
                dt = subsref(a,subsrefStruct);
                try
                    eval(['datetime(' rhs ', ''Format'', ''' strrep(dt.Format, '''', '''''') ''');']);
                    if row>size(a,1) || ~isequal(char(dt),rhsValue)
                        str =  [varname '.' colname '(' num2str(row) ') = ' rhs ';'];
                    end
                catch
                    error(message('MATLAB:datetime:InvalidFromVE'));
                end
            end
        else % e.g. x.Weight(10) = rhs;
            if isa(vardata,'categorical') 
                % Already enforced that rhsValue is a string for a categorical var
            elseif ~isscalar(rhsValue)
                if isempty(rhsValue)
                    error(message('MATLAB:datatypes:InvalidValue', colname));
                else
                    error(message('MATLAB:subsassignnumelmismatch'));
                end
            end
            subsrefStruct = struct('type',{'.','()'},'subs',{colname,{row}});
            % Generate code only if the new value is different
            if row>size(a,1) || ~isequal(subsref(a,subsrefStruct),rhsValue)
                % If the table has only one row and this edit expands the
                % number of rows, then specify the column as 1 to prevent
                % the table variable growing horizontally.
                if size(a,1)==1 && row>1 % force the var to grow vertically
                    str = [varname '.' colname '(' num2str(row) ',1) = ' rhs ';'];
                else 
                    str = [varname '.' colname '(' num2str(row) ') = ' rhs ';'];
                end
            end
        end
    end
elseif varIndices(end)==col % Editing one column beyond the table grid
    varName = matlab.internal.table.dfltVarNames(size(a,2)+1,true);
    varName = matlab.lang.makeUniqueStrings(varName,a.varnames,namelengthmax);
    
    % If rhs is a string or a numeric non-scalar, add the variable as a cell
    % array. Otherwise, add it as its native type.
    rhsValue = evalin('caller',rhs);
    if ~isempty(regexp(rhs,'^''.*''$','once')) || ischar(rhsValue) || ...
           iscell(rhsValue) || (isnumeric(rhsValue) && ~isscalar(rhsValue)) 
        str = [varname '.' varName  '{' num2str(row) ',1} = ' rhs ';'];
    else
        str = [varname '.' varName  '(' num2str(row) ',1) = ' rhs ';'];
    end
end


