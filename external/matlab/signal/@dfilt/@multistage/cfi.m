function fi = cfi(this)
%CFI   Return the information for the CFI.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

fi.Structure = get(this, 'FilterStructure');

try
    fi.Order = sprintf('%d', order(this));
catch
    fi.Order = 'Unknown';
end

fi.Stages    = sprintf('%d', nstages(this));

% [EOF]
