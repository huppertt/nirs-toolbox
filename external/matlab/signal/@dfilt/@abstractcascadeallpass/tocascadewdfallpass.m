function h = tocascadewdfallpass(this)
%TOCASCADEWDFALLPASS   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

c = coefficients(this);

h = dfilt.cascadewdfallpass(c{:});

% [EOF]
