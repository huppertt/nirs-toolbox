function [selected,locs] = getParamVal(inputs,kwds,argname,multok)
%getParamVal Select from a finite set of choices
%   [PVAL,PLOC] = getParamVal(PVAL,KWDS,PNAME) process input parameters
%   that can take values from a finite set. On input, PVAL is a candidate
%   parameter value (string), KWDS is a cell array of strings defining all
%   legal parameter values, and PNAME is the name of the parameter (used
%   only if needed in an error message). getParamVal tries to match PVAL
%   to the KWDS values, ignoring case and allowing unambigous abbreviations.
%   On output, PVAL is the matched value and PLOC is its index in KWDS.
%
%   [PVAL,PLOC] = getParamVal(PVAL,KWDS,PNAME,true) allows for multiple
%   selections. Here PVAL can be a cell array on input and is a cell
%   array on output.
%
%   Example:
%       onoff = getParamVal(pval,{'on' 'off'},'Display')


%   Copyright 2010-2011 The MathWorks, Inc.


if nargin<4
    multok = false;
end

if ~ischar(inputs) && ~iscellstr(inputs)
    if multok
        m = message('stats:internal:getParamVal:IllegalValues',argname);
        throwAsCaller(MException(m.Identifier,'%s',getString(m)));
    else
        m = message('stats:internal:getParamVal:IllegalValue',argname);
        throwAsCaller(MException(m.Identifier,'%s',getString(m)));
    end
end

% Convert to standard form for convenience
if ischar(inputs)
    inputs = cellstr(inputs);
else
    inputs = inputs(:);
end
if ischar(kwds)
    kwds = cellstr(kwds);
end

% Process each input value
n = length(inputs);
selected = cell(n,1);
locs = zeros(n,1);

% No empties, only one input allowed by default
if n==0
    m = message('stats:internal:getParamVal:EmptyValue',argname);
    throwAsCaller(MException(m.Identifier,'%s',getString(m)));
elseif n>1 && ~multok
    m = message('stats:internal:getParamVal:TooManyValues',argname);
    throwAsCaller(MException(m.Identifier,'%s',getString(m)));
end

% Process each input
for j=1:n
    txt = inputs{j};
    rows = find(strncmpi(txt,kwds,max(length(txt),1)));
    if isempty(rows)
        if length(kwds) < 6
            m = message('stats:internal:getParamVal:BadValueListChoices',txt,argname,liststrings(kwds));
        else % more than 5
            m = message('stats:internal:getParamVal:BadValue',txt,argname);
        end
        throwAsCaller(MException(m.Identifier,'%s',getString(m)));
    elseif ~isscalar(rows)
        k = strcmpi(txt,kwds(rows));
        if sum(k)==1
            rows = rows(k);
        else  % presumably this is empty unless the kwds input had repeats
            m = message('stats:internal:getParamVal:AmbiguousValue',txt,argname);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    end
    selected{j} = kwds{rows};
    locs(j) = rows;
end

% Return string unless multiple selections allowed
if ~multok
    selected = selected{1};
end

function s = liststrings(c)
n = numel(c);
fmt = '''%s''';
if n>1
    fmt = [fmt, repmat(', ''%s''',1,n-1)];
end
s = sprintf(fmt,c{:});
