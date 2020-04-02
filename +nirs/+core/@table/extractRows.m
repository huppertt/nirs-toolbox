function b = extractRows(t,rowIndices,extractCells)
%EXTRACTROW Retrieve one or more table rows as a 1-by- cell vector.

%   Copyright 2012 The MathWorks, Inc.

b = cell(1,t.nvars);
t_data = t.data;
for j = 1:t.nvars
    var_j = t_data{j};
    if ismatrix(var_j)
        b{j} = var_j(rowIndices,:); % without using reshape, may not have one
    else
        % Each var could have any number of dims, no way of knowing,
        % except how many rows they have.  So just treat them as 2D to get
        % the necessary rows, and then reshape to their original dims.
        sizeOut = size(var_j); sizeOut(1) = numel(rowIndices);
        b{j} = reshape(var_j(rowIndices,:), sizeOut);
    end
    if extractCells && iscell(b{j})
        b{j} = vertcat(b{j}{:});
    end
end
