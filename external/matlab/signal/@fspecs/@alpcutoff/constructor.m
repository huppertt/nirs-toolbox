function constructor(h,N,Wc)
%CONSTRUCTOR   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

if nargin > 1,
    h.FilterOrder = N;
end

if nargin > 2,
    h.Wcutoff = Wc;
end

% [EOF]
