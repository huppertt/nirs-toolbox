function [c,il,ir] = joinInnerOuter(a,b,leftOuter,rightOuter,leftKeyvals,rightKeyvals, ...
                                    leftVars,rightVars,leftVarnames,rightVarnames)
%JOININNEROUTER Common calculations for innerJoin and outerJoin.

% C is [A(IA,LEFTVARS) B(IB,RIGHTVARS)], where IA and IB are row indices into A
% and B computed for each row of C from LEFTKEYVALS and RIGHTKEYVALS.  These
% index vectors may include zeros indicating "no source row in A/B)" for some
% rows of C.

%   Copyright 2012 The MathWorks, Inc.

import matlab.internal.tableUtils.defaultarrayLike
import matlab.internal.tableUtils.repelem
import matlab.internal.tableUtils.coloncat

try
    % Sort each key.
    [lkeySorted,lkeySortOrd] = sort(leftKeyvals);
    [rkeySorted,rkeySortOrd] = sort(rightKeyvals);

    % Get unique key values and counts.   This also gives the beginning and end
    % of each block of key values in each.
    lbreaks = find(diff(lkeySorted));
    rbreaks = find(diff(rkeySorted));
    lstart = [ones(1,~isempty(leftKeyvals)); lbreaks+1]; % empty if left key is
    rstart = [ones(1,~isempty(rightKeyvals)); rbreaks+1]; % empty if right key is
    lend = [lbreaks; length(lkeySorted)];
    rend = [rbreaks; length(rkeySorted)];
    lunique = lkeySorted(lstart);
    runique = rkeySorted(rstart);
    luniqueCnt = lend - lstart + 1;
    runiqueCnt = rend - rstart + 1;
    clear lbreaks rbreaks lstart lend % clear some potentially large variables no longer needed

    % Use the "block nested loops" algorithm to determine how many times to
    % replicate each row of A and B.  Rows within each "constant" block of keys in
    % A will need to be replicated as many times as there are rows in the matching
    % block of B, and vice versa.  Rows of A that don't match anything in B, or
    % vice versa, get zero.  Rows of A will be replicated row-by-row; rows in B
    % will be replicated block-by-block.
    il = 1;
    ir = 1;
    leftElemReps = zeros(size(lunique));
    rightBlockReps = zeros(size(runique));
    while (il <= length(lunique)) && (ir <= length(runique))
        if lunique(il) < runique(ir)
            il = il + 1;
        elseif lunique(il) == runique(ir)
            leftElemReps(il) = runiqueCnt(ir);
            rightBlockReps(ir) = luniqueCnt(il);
            il = il + 1;
            ir = ir + 1;
        elseif lunique(il) > runique(ir)
            ir = ir + 1;
        else % one must have been NaN
            % NaNs get sorted to end; nothing else will match
            break;
        end
    end

    % Identify the rows of A required for an inner join: expand out the number of
    % replicates within each block to match against the (non-unique) sorted keys,
    % then replicate each row index the required number of times.
    leftElemReps = repelem(leftElemReps,luniqueCnt);
    il = repelem(1:length(lkeySorted),leftElemReps)';

    % Identify the rows of B required for an inner join: replicate the start and
    % end indices of each block of keys the required number of times, then create
    % a concatenation of those start:end expressions.
    rstart = repelem(rstart,rightBlockReps);
    rend = repelem(rend,rightBlockReps);
    ir = coloncat(rstart,rend)';
    clear rstart rend % clear some potentially large variables no longer needed

    % Translate back to the unsorted row indices.
    il = lkeySortOrd(il);
    ir = rkeySortOrd(ir);

    % If this is a left- or full-outer join, add the indices of the rows of A that
    % didn't match anything in B.  Add in zeros for the corresponding B indices.
    if leftOuter
        left = find(leftElemReps == 0);
        il = [il; lkeySortOrd(left)];
        ir = [ir; zeros(size(left))];
    end

    % If this is a right- or full-outer join, add the indices of the rows of B that
    % didn't match anything in A.  Add in zeros for the corresponding A indices.
    if rightOuter
        rightBlockReps = repelem(rightBlockReps,runiqueCnt);
        right = find(rightBlockReps == 0);
        il = [il; zeros(size(right))];
        ir = [ir; rkeySortOrd(right)];
    end

    % Now sort the whole thing by the key.  If this is an inner join, that's
    % already done.
    if leftOuter || rightOuter
        pos = (il > 0);
        Key = zeros(size(il));
        Key(pos) = leftKeyvals(il(pos)); % Rows that have an A key value
        Key(~pos) = rightKeyvals(ir(~pos)); % Rows with no A key value must have a B key
        [~,ord] = sort(Key);
        il = il(ord);
        ir = ir(ord);
    end

    % Create a new table by combining the specified variables from A and from B.
    c = cloneAsEmpty(a); % can't sensibly copy any Properties
    c.nrows = length(il);
    c.nvars = length(leftVars) + length(rightVars);
    c.varnames = [leftVarnames rightVarnames];
    c_data = cell(1,c.nvars);

    % Compute logical indices of where A'a and B's rows will go in C,
    % and the indices of which rows to pick out of A and B.
    ilDest = (il > 0); ilSrc = il(ilDest);
    irDest = (ir > 0); irSrc = ir(irDest);

    % Move data into C.
    a_data = a.data; c_nrows = c.nrows;
    for j = 1:length(leftVars)
        leftvar_j = a_data{leftVars(j)};
        szOut = size(leftvar_j); szOut(1) = c_nrows;
        cvar_j = defaultarrayLike(szOut,'Like',leftvar_j);
        cvar_j(ilDest,:) = leftvar_j(ilSrc,:);
        c_data{j} = reshape(cvar_j,szOut);
    end
    b_data = b.data;
    for j = 1:length(rightVars)
        rightvar_j = b_data{rightVars(j)};
        szOut = size(rightvar_j); szOut(1) = c_nrows;
        cvar_j = defaultarrayLike(szOut,'Like',rightvar_j);
        cvar_j(irDest,:) = rightvar_j(irSrc,:);
        c_data{length(leftVars) + j} = reshape(cvar_j,szOut);
    end
    c.data = c_data;

    % Carry over var descr and units from A and B, var by var.
    c.props.VariableDescriptions = catVarProps(a.props.VariableDescriptions,b.props.VariableDescriptions,leftVars,rightVars);
    c.props.VariableUnits = catVarProps(a.props.VariableUnits,b.props.VariableUnits,leftVars,rightVars);
catch ME
    throwAsCaller(ME)
end
