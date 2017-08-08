function [wmin, wmax] = getfreqinputs(hObj, unitcircle)
%GETFREQINPUTS Returns the frequency response inputs

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

error(nargchk(1,2,nargin,'struct'));
if nargin < 2
    unitcircle = 1;
end

% If getmaxfs returns [], one or more of the filters is normalized.  We
% then set the maxfs to 2 (since we normalize to pi, 2 represents 2pi).
wmax = getmaxfs(hObj);
if isempty(wmax)
    wmax = 2*pi;
end

switch unitcircle
    case 1,
        wmin = 0;
        wmax = wmax/2;
    case 2,
        wmin = 0;
    case 3,
        wmin = -wmax/2;
        wmax = wmax/2;
    otherwise
        wmin = [];
        wmax = [];
end

% [EOF]
