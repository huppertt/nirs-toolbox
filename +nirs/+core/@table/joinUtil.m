function [leftVars,rightVars,leftVarNames,rightVarNames,leftKeyVals,rightKeyVals,leftKeys,rightKeys] ...
                = joinUtil(a,b,type,leftTableName,rightTableName, ...
                           keys,leftKeys,rightKeys,leftVars,rightVars,keepOneCopy,supplied)
%JOINUTIL Common set-up for join, innerjoin, and outerjoin.

%   Copyright 2012 The MathWorks, Inc.
try
    rowNamesKey = false;
    if ~supplied.Keys
        if ~supplied.LeftKeys && ~supplied.RightKeys
            [leftKeys,rightKeys] = ismember(a.varnames,b.varnames);
            leftKeys = find(leftKeys);
            rightKeys = rightKeys(rightKeys>0);
            if isempty(leftKeys)
                error(message('MATLAB:table:join:CantInferKey'));
            end
        elseif ~supplied.LeftKeys || ~supplied.RightKeys
            error(message('MATLAB:table:join:MissingKeyVar'));
        else
            % Make sure the keys exist in both sides.
            leftKeys = getVarIndices(a,leftKeys);
            rightKeys = getVarIndices(b,rightKeys);
            if length(leftKeys) ~= length(rightKeys)
                error(message('MATLAB:table:join:UnequalNumKeyVars'));
            end
        end
    else
        if supplied.LeftKeys || supplied.RightKeys
            error(message('MATLAB:table:join:ConflictingInputs'));
        elseif isequal(keys,'RowNames') && strcmpi(type,'simple')
            % When using the row names no other keys are allowed,
            % and there'd be no point since the names are unique
            rowNamesKey = true;
            leftKeys = 0;
            rightKeys = 0;
        else
            leftKeys = getVarIndices(a,keys);
            rightKeys = getVarIndices(b,keys);
        end
    end

    % Use all vars from A and B by default, or use the specified vars.
    if supplied.LeftVariables
        leftVars = getVarIndices(a,leftVars);
        if length(unique(leftVars)) < length(leftVars)
            error(message('MATLAB:table:join:DuplicateVars'));
        end
    else
        leftVars = 1:a.nvars;
    end
    if supplied.RightVariables
        rightVars = getVarIndices(b,rightVars);
        if length(unique(rightVars)) < length(rightVars)
            error(message('MATLAB:table:join:DuplicateVars'));
        end
    else
        rightVars = 1:b.nvars;
        % Leave out B's keys for simple/inner joins, they are identical to A's keys.
        % Don't need to leave out keys if using row names.
        if ~rowNamesKey && (strcmpi(type,'simple') || strcmpi(type,'inner'))
            rightVars(rightKeys) = [];
        end
    end

    % Handle duplicate var names.
    leftVarNames = a.varnames(leftVars);
    rightVarNames = b.varnames(rightVars);
    [dups,ia,ib] = intersect(leftVarNames,rightVarNames);
    if supplied.KeepOneCopy && ~isempty(dups)
        [~,keepOneCopy] = intersect(dups,keepOneCopy);
        dropFromB = ib(keepOneCopy);
        rightVars(dropFromB) = [];
        rightVarNames(dropFromB) = [];
        [dups,ia,ib] = intersect(leftVarNames,rightVarNames);
    end
    if ~isempty(dups)
        % Uniqueify any duplicate var names.
        if isempty(leftTableName)
            leftTableName = getString(message('MATLAB:table:uistrings:JoinLeftVarSuffix'));
        end
        leftVarNames(ia) = strcat(leftVarNames(ia),['_' leftTableName]);
        if isempty(rightTableName)
            rightTableName = getString(message('MATLAB:table:uistrings:JoinRightVarSuffix'));
        end
        rightVarNames(ib) = strcat(rightVarNames(ib),['_' rightTableName]);
        % Don't allow the uniqueified names on either side to duplicate existing
        % names from either side
        vn = [leftVarNames rightVarNames];
        vn = matlab.lang.makeUniqueStrings(vn,ia,namelengthmax);
        vn = matlab.lang.makeUniqueStrings(vn,length(leftVarNames)+ib,namelengthmax);
        leftVarNames = vn(1:length(leftVarNames)); rightVarNames = vn(length(leftVarNames)+1:end);
    end

    if rowNamesKey
        [~,leftKeyVals,rightKeyVals] = intersect(a.rownames,b.rownames,'stable');
        if length(leftKeyVals) ~= a.nrows
            error(message('MATLAB:table:join:UnequalRowNames'));
        end

    else
        % Get the key var values, and check that they are scalar-valued or
        % vector-valued.
        leftKeyNames = a.varnames(leftKeys);
        rightKeyNames = b.varnames(rightKeys);
        leftKeyVals = a.data(leftKeys);
        rightKeyVals = b.data(rightKeys);
        if any(cellfun('ndims',leftKeyVals) > 2) || any(cellfun('ndims',rightKeyVals) > 2)
            error(message('MATLAB:table:join:NDKeyVar'));
        end

        % Convert possibly multiple keys to a single integer-valued key, taking on
        % comparable values across A and B.
        nkeys = length(leftKeys);
        leftlen = size(a,1);
        rightlen = size(b,1);
        lrkeys = zeros(leftlen+rightlen,nkeys);
        for j = 1:nkeys
            if size(leftKeyVals{j},2) ~= size(rightKeyVals{j},2) % already know these are 2-D
                 error(message('MATLAB:table:join:KeyVarSizeMismatch', leftKeyNames{j}, rightKeyNames{j}));
            elseif iscell(leftKeyVals{j}) ~= iscell(rightKeyVals{j})
                 error(message('MATLAB:table:join:KeyVarCellMismatch', leftKeyNames{j}, rightKeyNames{j}));
            end
            try
                lrkey_j = [leftKeyVals{j}; rightKeyVals{j}];
            catch me
                 error(message('MATLAB:table:join:KeyVarTypeMismatch', leftKeyNames{j}, rightKeyNames{j}));
            end
            if size(lrkey_j,2) > 1
                if isnumeric(lrkey_j) || islogical(lrkey_j) || ischar(lrkey_j)
                    [~,~,lrkeys(:,j)] = unique(lrkey_j,'rows');
                else
                    error(message('MATLAB:table:join:MulticolumnKeyVar', class( rightKeyVals )));
                end
            else
                try
                    [~,~,lrkeys(:,j)] = unique(lrkey_j);
                catch me
                    if strcmp(me.identifier,'MATLAB:UNIQUE:InputClass')
                        error(message('MATLAB:table:join:KeyVarNonStringError', leftKeyNames{j}, rightKeyNames{j}));
                    else
                        error(message('MATLAB:table:join:KeyVarUniqueError', leftKeyNames{j}, rightKeyNames{j}, me.message));
                    end
                end
            end
        end
        if nkeys > 1
            [~,~,lrkeys] = unique(lrkeys,'rows');
        end
        leftKeyVals = lrkeys(1:leftlen);
        rightKeyVals = lrkeys(leftlen+(1:rightlen));
    end
catch ME
    throwAsCaller(ME)
end