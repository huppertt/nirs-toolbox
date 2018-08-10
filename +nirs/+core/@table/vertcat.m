function t = vertcat(varargin)
%VERTCAT Vertical concatenation for tables.
%   T = VERTCAT(T1, T2, ...) vertically concatenates the tables T1, T2, ... .
%   Row names, when present, must be unique across tables.  VERTCAT fills
%   in default row names for the output when some of the inputs have names
%   and some do not.
%
%   Variable names for all tables must be identical except for order.  VERTCAT
%   concatenates by matching variable names.  VERTCAT assigns values for each
%   property (except for RowNames) in T using the first non-empty value for
%   the corresponding property in the arrays T1, T2, ... .
%
%   See also CAT, HORZCAT.

%   Copyright 2012-2014 The MathWorks, Inc.

try
    nrows = cellfun(@(t)size(t,1),varargin); % dispatch to overloaded size, not built-in
    t_nrows = sum(nrows);
    t_nvars = []; % don't know yet
    stillNeedVarNames = true;

    need.Descr = true; need.VarDescrs = true; need.VarUnits = true; need.UserData = true;

    b_is0x0 = false(nargin,1);
    for j = 1:nargin
        b = varargin{j};
        wasCell = iscell(b);
        if isequal(b,[]) % accept this as a valid "identity element"
            b = table;
        elseif wasCell
            b = cell2table(b); % default var names won't be used
        elseif ~isa(b,'nirs.core.table')
            error(message('MATLAB:table:vertcat:InvalidInput'));
        end

        if b.nvars==0 && b.nrows==0 % special cases to mimic built-in behavior
            b_is0x0(j) = true;
            % do nothing, b_data(j,:) already contains []'s
        elseif isempty(t_nvars) % first non-0x0 input
            t = b; % preserve the subclass
            t_nvars = t.nvars;
            [t_varnamesSorted,t_ord] = sort(t.varnames);
            t_reord(t_ord) = 1:t_nvars; %#ok<AGROW>
            t_rownames = t.rownames;
            stillNeedVarNames = wasCell; % use var names from a table, keep looking if a cell array

            % Set up row names if the first array has them
            if ~isempty(t_rownames)
                t_rownames = lengthenVar(t_rownames,t_nrows);
                dfltdRowNames = false(size(t_rownames));
            end

            b_data = cell(nargin,t_nvars);
            b_data(j,:) = b.data;
        elseif b.nvars ~= t_nvars
            error(message('MATLAB:table:vertcat:SizeMismatch'));
        else
            if wasCell
                % Assign positionally
                b_data(j,:) = b.data;
            else % was always a table
                if stillNeedVarNames
                    t.varnames = b.varnames;
                    [t_varnamesSorted,t_ord] = sort(t.varnames);
                    t_reord(t_ord) = 1:t_nvars; %#ok<AGROW>
                    stillNeedVarNames = false;
                    b_data(j,:) = b.data;
                else
                    %[tf,b_reord] = ismember(t.varnames,b.varnames);
                    [b_varnamesSorted,b_ord] = sort(b.varnames);
                    tf = strcmp(t_varnamesSorted,b_varnamesSorted);
                    if ~all(tf)
                        error(message('MATLAB:table:vertcat:UnequalVarNames'));
                    end
                    b_data(j,:) = b.data(b_ord(t_reord));
                end
            end

            b_rownames = b.rownames;
            rows_j = sum(nrows(1:j-1)) + (1:nrows(j));
            if ~isempty(t_rownames) && ~isempty(b_rownames)
                % Save row names from this array
                checkDuplicateNames(b_rownames,t_rownames,'rownames');
                t_rownames(rows_j) = b_rownames;
            elseif ~isempty(t_rownames) % && isempty(b_rownames)
                % Create default row names for this array
                t_rownames(rows_j) = matlab.internal.table.dfltRowNames(rows_j);
                dfltdRowNames(rows_j) = true;
            elseif ~isempty(b_rownames) % && isempty(t_rownames)
                % Set up default row names for the previous arrays if we first
                % encounter row names part-way through
                t_rownames = cell(t_nrows,1);
                dfltdRowNames = false(size(t_rownames));
                rows_0 = 1:(rows_j(1)-1);
                t_rownames(rows_0) = matlab.internal.table.dfltRowNames(rows_0);
                dfltdRowNames(rows_0) = true;
                t_rownames(rows_j) = b_rownames;
            end
        end

        % Find the first non-empty property values.
        b_props = b.props;
        if need.Descr && ~isempty(b_props.Description)
            t.props.Description = b_props.Description;
            need.Descr = false;
        end
        if need.VarDescrs && ~isempty(b_props.VariableDescriptions)
            t.props.VariableDescriptions = b_props.VariableDescriptions;
            need.VarDescrs = false;
        end
        if need.VarUnits && ~isempty(b_props.VariableUnits)
            t.props.VariableUnits = b_props.VariableUnits;
            need.VarUnits = false;
        end
        if need.UserData && ~isempty(b_props.UserData)
            t.props.UserData = b_props.UserData;
            need.UserData = false;
        end
        % DimensionNames is always non-empty, will use a's.
    end

    if all(b_is0x0) % all inputs were empty
        t = varargin{1}; % use the first input
        t_nvars = 0;
        t_rownames = {};
    end

    t_data = cell(1,t_nvars);
    for i = 1:t_nvars
        try
            t_data{i} = vertcat(b_data{:,i}); % []'s are dropped
        catch ME
            m = message('MATLAB:table:vertcat:VertcatMethodFailed',t.varnames{i});
            throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
        end
        % Something went badly wrong with whatever vertcat method was called.
        if size(t_data{i},1) ~= t_nrows
            % One reason for this is concatenation of a cell variable with a non-cell
            % variable, which adds only a single cell to the former, containing the
            % latter.  Check for cell/non-cell only after calling vertcat to allow
            % overloads such as categorical that can vertcat cell/non-cell sensibly.
            cells = cellfun('isclass',b_data(~b_is0x0,i),'cell');
            if any(cells) && ~all(cells)
                error(message('MATLAB:table:vertcat:VertcatCellAndNonCell', t.varnames{i}));
            else
                error(message('MATLAB:table:vertcat:VertcatWrongLength', t.varnames{i}));
            end
        end
    end
    t.data = t_data;
    t.nrows = t_nrows;
    t.nvars = t_nvars;
    if ~isempty(t_rownames)
        % make sure default names don't conflict with supplied names
        t.rownames = matlab.lang.makeUniqueStrings(t_rownames,dfltdRowNames,namelengthmax);
    end
catch ME
    throwAsCaller(ME)
end
