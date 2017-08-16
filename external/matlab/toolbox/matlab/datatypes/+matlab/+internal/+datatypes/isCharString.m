function tf = isCharString(s,allowEmpty)
%ISCHARSTRING True for a string
%   T = ISCHARSTRING(S) returns true if S is a 1xN character vector
%   for N >= 0, or the 0x0 char array ''.
%
%   T = ISCHARSTRING(S,FALSE) returns true only if S is a non-empty
%   character vector.

%   Copyright 2012-2016 The MathWorks, Inc.

if nargin < 2 || allowEmpty
    tf = ischar(s) && (isrow(s) || isequal(s,''));
else
    tf = ischar(s) && isrow(s) && ~all(isspace(s));
end
