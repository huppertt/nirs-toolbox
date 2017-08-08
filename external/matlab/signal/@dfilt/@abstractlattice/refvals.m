function rcvals = refvals(this)
%REFVALS   Reference coefficient values.
%This should be a private method.
%   The values are returned in a cell array.

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.

rcnames = refcoefficientnames(this);

rcvals = get(this,rcnames);

% [EOF]