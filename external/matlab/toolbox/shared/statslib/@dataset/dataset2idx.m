function [ainds,binds] = dataset2idx(a,avars,b,bvars)
% DATASET2IDX Create index vectors from dataset variables.

%   Copyright 2012 The MathWorks, Inc.


paired = (nargin == 4);

nvars = length(avars); % == length(bvars) for paired
a_data = a.data;
anobs = a.nobs;
ainds = zeros(anobs,nvars);

if paired
    if ~isequal(b.varnames(bvars),a.varnames(avars))
        error(message('stats:dataset:setmembership:DifferentVars'));
    end
    bnobs = b.nobs;
    binds = zeros(bnobs,nvars);
    b_data = b.data;
end

for j = 1:nvars
    if paired
        avar_j = a_data{avars(j)}; avar_j = avar_j(:,:);
        bvar_j = b_data{bvars(j)}; bvar_j = bvar_j(:,:);
        try
            var_j = [avar_j; bvar_j];
        catch ME
            m = message('stats:dataset:setmembership:VarVertcatMethodFailed',a.varnames{avars(j)});
            throwAsCaller(addCause(MException(m.Identifier,'%s',getString(m)), ME));
        end
    else
        var_j = a_data{avars(j)}; var_j = var_j(:,:);
    end
    
    try
        % Use 'rows' for this variable's union method if the var is
        % not a single column.
        if iscolumn(var_j)
            [~,~,inds_j] = unique(var_j,'sorted');
        else
            [~,~,inds_j] = unique(var_j,'sorted','rows');
        end
    catch ME
        m = message('stats:dataset:setmembership:VarUniqueMethodFailed',a.varnames{avars(j)});
        throwAsCaller(addCause(MException(m.Identifier,'%s',getString(m)), ME));
    end
    if length(inds_j) ~= length(var_j)
        m = message('stats:dataset:setmembership:VarUniqueMethodFailedNumRows', a.varnames{avars(j)});
        throwAsCaller(MException(m.Identifier,'%s',getString(m)));
    end
    ainds(:,j) = inds_j(1:anobs,1);
    if paired, binds(:,j) = inds_j((anobs+1):end,1); end
end
