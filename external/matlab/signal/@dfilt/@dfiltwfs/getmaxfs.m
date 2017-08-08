function [fs, xunits] = getmaxfs(h)
%MAXFS Returns the maximum fs

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

fs = get(h, 'Fs');

if iscell(fs),
    fs = max([fs{:}]);
end

if nargout > 1,
    if isempty(fs),
        xunits = 'rad/sample';
    else
        [fs, m, xunits] = engunits(fs);
        xunits          = [xunits 'Hz'];
    end
end

% [EOF]
