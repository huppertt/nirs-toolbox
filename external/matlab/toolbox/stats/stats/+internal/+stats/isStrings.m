function [tf,asCellStr] = isStrings(s,requireCell)
% ISSTRINGS Require a char row vector, or '', or a cell array of same
%   T = ISSTRINGS(S) returns true if S is a string, as defined by ISSTRING, or
%   if S is a cell array containing strings, and false otherwise.
%
%   T = ISSTRINGS(S,REQUIRECELL), when REQUIRECELL is true, returns true only
%   when S is a cell array containing strings, as defined by ISSTRING, and
%   false otherwise.
%
%   When S is a single char string, [T,ASCELLSTR] = ISSTRINGS(S) also returns
%   {S}.  If S is a cell array containing strings, ISSTRINGS returns S itself
%   in ASCELLSTR.  Otherwise, ISSTRINGS returns {} in ASCELLSTR.
%   
%   See also ISSTRING, STRINGS.


%   Copyright 2011 The MathWorks, Inc.

stringTest = @(s) ischar(s) && ( (isvector(s) && (size(s,1) == 1)) || all(size(s)==0) );
if iscell(s)
    tf = all(cellfun(stringTest,s,'UniformOutput',true));
    if nargout > 1
        if tf
            asCellStr = s;
        else
            asCellStr = {};
        end
    end
elseif nargin < 2 || ~requireCell
    tf = stringTest(s);
    if nargout > 1
        if tf
            asCellStr = {s};
        else
            asCellStr = {};
        end
    end
else
    tf = false;
    asCellStr = {};
end
end
