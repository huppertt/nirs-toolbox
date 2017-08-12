function [tf,asChar] = isString(s,allowCell)
% ISSTRING Require a char row vector, or ''
%   T = ISSTRING(S) returns true if S is a 1-by-M character vector, where M
%   can be zero, or if S is the empty string '', and returns false otherwise.
%   
%   T = ISSTRING(S,ALLOWCELL), when ALLOWCELL is true, returns true if S is
%   either a string or a scalar cell array containing a single string, and
%   returns false otherwise.
%   
%   When T is true, [T,ASCHAR] = ISSTRING(S,ALLOWCELL) also returns either S,
%   or S{1} if S is a scalar cell containing a single string.  When T is
%   false, ISSTRING returns '' in ASCELLSTR.
%   
%   See also STRINGS.


%   Copyright 2011 The MathWorks, Inc.

if nargin > 1 && allowCell && isscalar(s) && iscell(s)
    s = s{1};
end

tf = ischar(s) && ...  
        (   (isvector(s) && (size(s,1) == 1)) ...  % 'typical string'
         || all(size(s)==0) );                     % ''
if nargout > 1
    if tf
        asChar = s;
    else
        asChar = '';
    end
end
