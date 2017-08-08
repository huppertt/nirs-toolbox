function [dindx, qindx] = getfiltindx(h)
%GETFILTINDX Returns the indices of the dfilts

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

G = get(h, 'Filter');
if ~iscell(G), G = {G}; end

qindx     = [];
dindx     = [];
otherindx = [];

for n = 1:length(G),
    if isquantized(G{n})
        qindx = [qindx n];
    else
        dindx = [dindx n];
    end
end

% [EOF]
