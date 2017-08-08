function fi = cfi(this)
%CFI   Return the information for 

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

fi.Structure = get(this, 'FilterStructure');
fi.Order     = sprintf('%d', order(this));
fi.Sections  = sprintf('%d', nsections(this));

% [EOF]
