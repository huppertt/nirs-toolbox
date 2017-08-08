function snewstr = disp4filtstates(Hd,s)
%DISP4FILTSTATES Utility for objects which use a FILTSTATES.DFIIR object.
  
% This should be a private method.

%   Author: P. Costa
%   Copyright 1988-2004 The MathWorks, Inc.

% Code to display the numerator and denominator states
snewstr = evalc('s');

% Offset required to align the N in Numerator with the D in Denominator in
% the DFILT display
x = 22;
snewstr = strrep(snewstr, '[1x1 filtstates.dfiir]', dispstr(Hd.States, x));

% Remove the first line since it reads 'snew = \n'
sndx = strfind(snewstr, 's =');

if strcmpi(get(0, 'formatspacing'), 'loose')
    snewstr(1:sndx+4) = [];
else
    snewstr(1:sndx+3) = [];
end

% Look for extra new line feeds, we want to get rid of them all. Using
% regexp so that the cases when format compact or loose are used, the
% display continues to work.
sndx = min(regexp(snewstr, '\w'));  % First occurrence of a word
sndx = max(find(snewstr(1:sndx-1) == char(10)));
snewstr(1:sndx) = [];

% [EOF]
