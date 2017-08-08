function s = isfssame(hObj)
%ISFSSAME Returns true if the fs is the same for all the filters

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

allfs    = get(hObj, 'Fs');
if ~iscell(allfs),
    s = true;
else
    allfs = [allfs{:}];
    if isempty(allfs),
        s = true;
    elseif any(diff(allfs))
        s = false;
    else
        s = true;
    end
end

% [EOF]
