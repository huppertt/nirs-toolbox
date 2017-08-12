function a = horzcat(varargin)
%HORZCAT Horizontal concatenation for dataset arrays.
%   DS = HORZCAT(DS1, DS2, ...) horizontally concatenates the dataset arrays
%   DS1, DS2, ... .  You may concatenate dataset arrays that have duplicate
%   variable names, however, the variables must contain identical data, and
%   HORZCAT includes only one copy of the variable in the output dataset.
%
%   Observation names for all dataset arrays that have them must be identical
%   except for order. HORZCAT concatenates by matching observation names when
%   present, or by position for datasets that do not have observation names.
%
%   See also DATASET/CAT, DATASET/VERTCAT.

%   Copyright 2006-2011 The MathWorks, Inc.


b = varargin{1};
if isequal(b,[]) % accept this as a valid "identity element"
    b = dataset;
elseif ~isa(b,'dataset')
    error(message('stats:dataset:horzcat:InvalidInput'));
end
a = b;
if ~isempty(a.obsnames)
    [a_sortobsnames,a_obsord] = sort(a.obsnames);
end
for j = 2:nargin
    b = varargin{j};
    if isequal(b,[]) % accept this as a valid "identity element"
        b = dataset;
    elseif ~isa(b,'dataset')
        error(message('stats:dataset:horzcat:InvalidInput'));
    end
    
    % some special cases to mimic built-in behavior
    if a.nvars==0 && a.nobs==0
        a = b;
        if ~isempty(a.obsnames)
            [a_sortobsnames,a_obsord] = sort(a.obsnames);
        end
        continue;
    elseif b.nvars==0 && b.nobs==0
        % do nothing
        continue;
    elseif a.nobs ~= b.nobs
        error(message('stats:dataset:horzcat:SizeMismatch'));
    end
    
    [dups,dupBVars] = checkduplicatenames(b.varnames,a.varnames,'varnames');
    if dups
        for i = find(dupBVars)
            if ~isequal(a.data{getvarindices(a,b.varnames(i))},b.data{i})
                dup = b.varnames{i};
                error(message('stats:dataset:horzcat:DuplicateVarnames',dup));
            end
        end
        uniqBVars = find(~dupBVars);
    else
        uniqBVars = 1:b.nvars;
    end
    
    if ~isempty(a.obsnames) && ~isempty(b.obsnames)
        [b_sortobsnames,b_obsord] = sort(b.obsnames);
        if ~all(strcmp(a_sortobsnames,b_sortobsnames))
            error(message('stats:dataset:horzcat:UnequalObsnames'));
        end
        b_reord(a_obsord) = b_obsord; %#ok<AGROW>, full reassignment each time
        a.data = horzcat(a.data, cell(1,length(uniqBVars)));
        for i = 1:length(uniqBVars)
            bVar = b.data{uniqBVars(i)};
            sizeOut = size(bVar);
            a.data{a.nvars+i} = reshape(bVar(b_reord,:),sizeOut);
        end
    else
        if isempty(a.obsnames) && ~isempty(b.obsnames)
            a.obsnames = b.obsnames;
            [a_sortobsnames,a_obsord] = sort(a.obsnames);
        end
        a.data = horzcat(a.data, b.data(uniqBVars));
    end
    % We've already weeded out duplicate var names, either by erroring or by
    % ignoring the second copy, so there's no need to uniqueify var names.
    AVars = 1:a.nvars;
    a.nvars = a.nvars + length(uniqBVars);
    a.varnames = horzcat(a.varnames, b.varnames(uniqBVars));
    
    % Concatenate var descriptions and units.
    a.props.VarDescription = catVarProps(a.props.VarDescription,b.props.VarDescription,AVars,uniqBVars);
    a.props.Units = catVarProps(a.props.Units,b.props.Units,AVars,uniqBVars);
end
