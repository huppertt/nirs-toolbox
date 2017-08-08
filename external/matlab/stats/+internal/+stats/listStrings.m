function s = listStrings(c)
% LISTSTRINGS Create a list from a cell array of strings.
%    T = LISTSTRINGS(C) takes the strings in the cell array C and creates a
%    comma-separated list from them.


%   Copyright 2011 The MathWorks, Inc.

s = sprintf(', %s',c{:});
s(1:2) = [];  % remove leading comma
