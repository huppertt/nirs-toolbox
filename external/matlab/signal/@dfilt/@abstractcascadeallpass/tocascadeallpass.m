function h = tocascadeallpass(this)
%TOCASCADEALLPASS   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

c = coefficients(this);

h = dfilt.cascadeallpass(c{:});

% [EOF]
