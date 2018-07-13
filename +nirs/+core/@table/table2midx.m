function [ainds,binds] = table2midx(a,avars,b,bvars)
% TABLE2MIDX Create multi-index matrices from table variables.

%   Copyright 2012-2014 The MathWorks, Inc.

import matlab.internal.tableUtils.matricize

try
    paired = (nargin == 4);

    nvars = length(avars); % == length(bvars) for paired
    a_data = a.data;
    anrows = a.nrows;
    ainds = zeros(anrows,nvars);

    if paired
        if ~isequal(b.varnames(bvars),a.varnames(avars))
            error(message('MATLAB:table:setmembership:DifferentVars'));
        end
        bnrows = b.nrows;
        binds = zeros(bnrows,nvars);
        b_data = b.data;
    end

    for j = 1:nvars
        if paired
            avar_j = matricize(a_data{avars(j)});
            bvar_j = matricize(b_data{bvars(j)});
            try
                % Call union function/method first to enforce rules about mixed types
                % that vertcat is not strict about, e.g. single and int32
                union(avar_j([],:),bvar_j([],:));
                var_j = [avar_j; bvar_j];
            catch ME
                m = message('MATLAB:table:setmembership:VarVertcatMethodFailed',a.varnames{avars(j)});
                throwAsCaller(addCause(MException(m.Identifier,'%s',getString(m)), ME));
            end
        else
            var_j = matricize(a_data{avars(j)});
        end

        % unique won't work right on multi-column cellstrs catch these here to avoid
        % the 'rows' warning which would be followed by an error
        if iscell(var_j) && ~iscolumn(var_j)
            m = message('MATLAB:table:VarUniqueMethodFailedCellRows',a.varnames{avars(j)});
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end

        try
            % Use 'rows' for this variable's unique method if the var is
            % not a single column. Multi-column cellstrs already weeded out.
            %
            % For categorical variables, the indices created here _will_ account for
            % categories that are not actually present in the data -- the indices
            % should not be assumed to be contiguous.
            if iscolumn(var_j)
                [~,~,inds_j] = unique(var_j,'sorted');
            else
                [~,~,inds_j] = unique(var_j,'sorted','rows');
            end
        catch ME
            m = message('MATLAB:table:VarUniqueMethodFailed',a.varnames{avars(j)});
            throwAsCaller(addCause(MException(m.Identifier,'%s',getString(m)), ME));
        end
        if length(inds_j) ~= size(var_j,1)
            m = message('MATLAB:table:VarUniqueMethodFailedNumRows',a.varnames{avars(j)});
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
        ainds(:,j) = inds_j(1:anrows,1);
        if paired, binds(:,j) = inds_j((anrows+1):end,1); end
    end
catch ME
    throwAsCaller(ME)
end
