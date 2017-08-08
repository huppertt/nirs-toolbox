function S = quantizestates(q,S)
%QUANTIZESTATES   

%   Author(s): V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

if strcmpi(class(S),'filtstates.dfiir'),
    S.Numerator = double(S.Numerator);
    S.Denominator = double(S.Denominator);
else
    S = double(S);
end

% [EOF]
