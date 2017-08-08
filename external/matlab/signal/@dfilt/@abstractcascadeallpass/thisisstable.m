function b = thisisstable(this)
%THISISSTABLE   True if the object is stable.

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

Hd = dispatch(this);

b = isstable(Hd);


% [EOF]
