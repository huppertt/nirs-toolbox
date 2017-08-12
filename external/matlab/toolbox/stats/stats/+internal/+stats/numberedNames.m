function names = numberedNames(baseName,n)
% NUMBEREDNAMES Array of names formed by appending integer to prefix.
%    NAMES = NUMBEREDNAMES(PREFIX,N) accepts a string PREFIX and a
%    vector N of integers, and returns a cell array NAMES containing
%    strings formed by concatenating the string to the integers
%    from in N. If N is a vector of integers, NUMBEREDNAMES uses
%    the integers from the vector N.


%   Copyright 2011 The MathWorks, Inc.

if ~internal.stats.isString(baseName)
    error(message('stats:internal:numberedNames:InvalidPrefix'));
end

if isempty(n)
    names = cell(1,0);
else
    if ~isvector(n) || ~internal.stats.isIntegerVals(n,0)
        error(message('stats:internal:numberedNames:InvalidN'));
    end
    names = strcat({baseName},num2str(n(:),'%-d'))';
end
