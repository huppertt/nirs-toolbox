function t = strCollapse(s,sep)
%STRCOLLAPSE Concatenate strings into a single string.
%   T = STRCOLLAPSE(S) concatenates all the strings in S into a single string.
%   S is either a character array, where each row is a string, or a cell array
%   of strings.  T is a character row vector.
%
%   T = STRCOLLAPSE(S,SEP) concatenates the strings in S into a single string,
%   separated by the string SEP.
%
%   STRCOLLAPSE ignores trailing ASCII white space characters and omits all 
%   such characters from the output.  White space characters in ASCII are 
%   space, newline, carriage return, tab, vertical tab, or form-feed 
%   characters, all of which return a TRUE response from the MATLAB ISSPACE
%   function.
%
%   Example
%       strcollapse({'Red','Yellow','Green','Blue'})
%   returns
%       'RedYellowGreenBlue'
%   while
%       strcollapse({'Red','Yellow','Green','Blue'},',')
%   returns
%       'Red,Yellow,Green,Blue'
%
%   See also STRCAT, CELLSTR.


%   Copyright 2011 The MathWorks, Inc.

if ischar(s)
    s = cellstr(s);
end
if ~internal.stats.isStrings(s)
    error(message('stats:internal:strCollapse:InvalidStrings'));
end

if nargin > 1
    if ~internal.stats.isString(sep)
        error(message('stats:internal:strCollapse:InvalidSeparator'));
    end
    s = [s(:)'; repmat({sep},1,numel(s))];
    if ~isempty(s)
        s{end} = ''; 
    end
end
if isempty(s)
    t = '';
else
    t = cat(2,s{:}); % strcat would remove trailing spaces
end
