function [str,msg] = variableEditorSetDataCode(a,varname,row,col,rhs)
%   This function is undocumented and will change in a future release

% Generate MATLAB command to edit the content of a cell to the specified
% rhs.

%   Copyright 2011-2013 The MathWorks, Inc.

msg = '';
str = '';

[colNames,varIndices,colClasses] = variableEditorColumnNames(a);

if col>varIndices(end)
    msg = getString(message('stats:dataset:VarEditorIndexOverlow'));
    col = varIndices(end);
end

% Find the left hand column of the dataset variable that contains col
j = find(varIndices>=col+1,1,'first')-1;

if ~isempty(j)
    colname = colNames{j};
    vardata = a.data{j};
    isGroupColumn = varIndices(j+1)>varIndices(j)+1;
    isCatColumn = isa(vardata,'categorical');
    if isCatColumn
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
        else
            rhs = '[]';
        end
    end
    rhsValue = evalin('caller',rhs);
    
    % Generate code only when rhs is a new value. This avoids cluttering the
    % command window when the dataset element would not change. Upstream code
    % catches a rhs that is literally the same, here we catch cases where it is
    % an expression that evaluates to the same value.
        
    % If this is a grouped column, an additonal numeric sub-index is needed
    % to define which column within the dataset is being edited. A char matrix
    % is treated as one column, not as grouped columns.
    if isGroupColumn
        relcol = col-varIndices(j)+1;
        if iscell(vardata) % e.g. x.BloodPressure{10,2} = rhs;
            % Generate code only if the new value is different
            subsrefStruct = struct('type',{'.','{}'},'subs',{colname,{row,relcol}});
            if row>size(a,1) || ~isequal(subsref(a,subsrefStruct),rhsValue)
                str =  [varname '.' colname '{' num2str(row) ', ' num2str(relcol) '} = ' rhs ';'];
            end
        else % e.g. x.BloodPressure(10,2) = rhs;
            if isa(vardata,'categorical') 
               % Already enforced that rhsValue is a string for a categorical var
            elseif ~isscalar(rhsValue)
                error(message('MATLAB:subsassignnumelmismatch'));
            end
            % Generate code only if the new value is different
            subsrefStruct = struct('type',{'.','()'},'subs',{colname,{row,relcol}});
            if row>size(a,1) || ~isequal(subsref(a,subsrefStruct),rhsValue)
                str =  [varname '.' colname '(' num2str(row) ','  num2str(relcol) ') = ' rhs ';'];
            end
        end
    else
        if iscell(vardata) % e.g. x.Weight{10} = rhs;
            % Generate code only if the new value is different
            subsrefStruct = struct('type',{'.','{}'},'subs',{colname,{row}});
            if row>size(a,1) || ~isequal(subsref(a,subsrefStruct),rhsValue)
                if size(a,1)==1 && row>1 % force the var to grow vertically
                    str =  [varname '.' colname '{' num2str(row) ',1} = ' rhs ';'];
                else
                    str =  [varname '.' colname '{' num2str(row) '} = ' rhs ';'];
                end
            end
        elseif strcmp('char',colClasses{j})
            % Generate code only if the new value is different
            subsrefStruct = struct('type',{'.','()'},'subs',{colname,{row,':'}});
            if row>size(a,1) || ~isequal(subsref(a,subsrefStruct),rhsValue)
                str = [varname '.' colname '(' num2str(row) ',:) = ' rhs ';'];
            end
        else % e.g. x.Weight(10) = rhs;
            if isa(vardata,'categorical') 
               % Already enforced that rhsValue is a string for a categorical var
            elseif ~isscalar(rhsValue)
                error(message('MATLAB:subsassignnumelmismatch'));
            end
            subsrefStruct = struct('type',{'.','()'},'subs',{colname,{row}});
            % Generate code only if the new value is different
            if row>size(a,1) || ~isequal(subsref(a,subsrefStruct),rhsValue)
                % If the dataset has only one row and this edit expands the
                % number of rows, then specify the column as 1 to prevent
                % the dataset variable growing horizontally.
                if size(a,1)==1 && row>1 % force the var to grow vertically
                    str =  [varname '.' colname '(' num2str(row) ',1) = ' rhs ';'];
                else 
                    str =  [varname '.' colname '(' num2str(row) ') = ' rhs ';'];
                end
            end
        end
    end
elseif varIndices(end)==col % Editing one column beyond the table grid
    varName = dfltvarnames(size(a,2)+1,true);
    varName = matlab.lang.makeUniqueStrings(varName,a.varnames,namelengthmax);
    
    % If rhs is a string or a numeric non-scalar, add the variable as a cell
    % array. Otherwise, add it as its native type.
    rhsValue = evalin('caller',rhs);
    if ~isempty(regexp(rhs,'^''.*''$','once')) || ischar(rhsValue) || ...
           iscell(rhsValue) || (isnumeric(rhsValue) && ~isscalar(rhsValue)) 
        str =  [varname '.' varName  '{' num2str(row) ',1} = ' rhs ';'];
    else
        str =  [varname '.' varName  '(' num2str(row) ',1) = ' rhs ';'];
    end
end


