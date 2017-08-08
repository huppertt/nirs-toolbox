function S = validatestatesobj(q, S)
%VALIDATESTATESOBJ   

%   Author(s): V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

S.Numerator = validatestates(q,S.Numerator);
S.Denominator = validatestates(q,S.Denominator);

% [EOF]
