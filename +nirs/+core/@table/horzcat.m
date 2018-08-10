function t = horzcat(varargin)
%HORZCAT Horizontal concatenation for tables.
%   T = HORZCAT(T1, T2, ...) horizontally concatenates the tables T1, T2,
%   ... .  All inputs must have unique variable names.
%
%   Row names for all tables that have them must be identical except for order.
%   HORZCAT concatenates by matching row names when present, or by position for
%   tables that do not have row names.  HORZCAT assigns values for the
%   Description and UserData properties in T using the first non-empty value
%   for the corresponding property in the arrays T1, T2, ... .
%
%   See also CAT, VERTCAT, JOIN.

%   Copyright 2012-2014 The MathWorks, Inc.
try
    dfltdVarNames = false(1,0);

    b = varargin{1};
    if isequal(b,[]) % accept this as a valid "identity element"
        b = table;
    elseif iscell(b)
        jvars = 1:size(b,2);
        dfltdVarNames(jvars) = true;
        b = cell2table(b,'VariableNames',matlab.internal.table.dfltVarNames(jvars));
    elseif ~isa(b,'nirs.core.table')
        error(message('MATLAB:table:horzcat:InvalidInput'));
    end
    t = b;
    t_data = t.data;
    t_nvars = t.nvars;
    t_nrows = t.nrows;
    t_varnames = t.varnames;
    t_rownames = t.rownames;
    t_props = t.props;
    if ~isempty(t.rownames)
        [t_sortrownames,t_rowOrder] = sort(t.rownames);
    end

    need.Descr = isempty(t_props.Description); need.UserData = isempty(t_props.UserData);

    for j = 2:nargin
        b = varargin{j};
        wasCell = iscell(b);
        if isequal(b,[]) % accept this as a valid "identity element"
            b = table;
        elseif wasCell
            jvars = t_nvars+(1:size(b,2));
            dfltdVarNames(jvars) = true;
            b = cell2table(b,'VariableNames',matlab.internal.table.dfltVarNames(jvars));
        elseif ~isa(b,'nirs.core.table')
            error(message('MATLAB:table:horzcat:InvalidInput'));
        end

        % some special cases to mimic built-in behavior
        if t_nvars==0 && t.nrows==0
            t = b;
            t_data = t.data;
            t_nvars = t.nvars;
            t_nrows = t.nrows;
            t_varnames = t.varnames;
            t_rownames = t.rownames;
            t_props = t.props;
            need.Descr = isempty(t_props.Description); need.UserData = isempty(t_props.UserData);
            if ~isempty(t_rownames)
                [t_sortrownames,t_rowOrder] = sort(t_rownames);
            end
            continue;
        elseif b.nvars==0 && b.nrows==0
            % do nothing
            continue;
        elseif t_nrows ~= b.nrows
            error(message('MATLAB:table:horzcat:SizeMismatch'));
        end

        if ~isempty(t_rownames) && ~isempty(b.rownames)
            [b_sortrownames,b_rowOrder] = sort(b.rownames);
            if ~all(strcmp(t_sortrownames,b_sortrownames))
                error(message('MATLAB:table:horzcat:UnequalRowNames'));
            end
            b_reord(t_rowOrder) = b_rowOrder; %#ok<AGROW>, full reassignment each time
            t_data = horzcat(t_data, cell(1,b.nvars)); %#ok<AGROW>
            for i = 1:b.nvars
                bVar = b.data{i};
                sizeOut = size(bVar);
                t_data{t_nvars+i} = reshape(bVar(b_reord,:),sizeOut);
            end
        else
            if isempty(t_rownames) && ~isempty(b.rownames)
                t_rownames = b.rownames;
                [t_sortrownames,t_rowOrder] = sort(t_rownames);
            end
            t_data = horzcat(t_data, b.data); %#ok<AGROW>
        end

        % Concatenate var descriptions and units.
        AVars = 1:t_nvars;
        BVars = 1:b.nvars;
        b_props = b.props;
        t_props.VariableDescriptions = catVarProps(t_props.VariableDescriptions,b_props.VariableDescriptions,AVars,BVars);
        t_props.VariableUnits = catVarProps(t_props.VariableUnits,b_props.VariableUnits,AVars,BVars);

        % Find the first non-empty property values for the remaining properties.
        if need.Descr && ~isempty(b_props.Description)
            t_props.Description = b_props.Description;
            need.Descr = false;
        end
        if need.UserData && ~isempty(b_props.UserData)
            t_props.UserData = b_props.UserData;
            need.UserData = false;
        end

        % Will still need to check for duplicate var names, do that afterwards.
        % Don't try to recognize and take one copy of variables that have the same
        % name and data, duplicate variable names is an error regardless.
        t_nvars = t_nvars + b.nvars;
        t_varnames = horzcat(t_varnames, b.varnames); %#ok<AGROW>
    end
    t.data = t_data;
    t.nvars = t_nvars;
    % t.nrows stays the same

    % Make sure that variable names created automatically for cell arrays don't
    % conflict with variable names from tables.  Then error if any variables have
    % the same name.
    if any(dfltdVarNames)
        if length(dfltdVarNames) < length(t_varnames)
            dfltdVarNames(length(t_varnames)) = false;
        end
        t_varnames = matlab.lang.makeUniqueStrings(t_varnames,dfltdVarNames,namelengthmax);
    end
    checkDuplicateNames(t_varnames,'varnames');

    t.varnames = t_varnames;
    t.rownames = t_rownames;
    t.props = t_props;
catch ME
    throwAsCaller(ME)
end
