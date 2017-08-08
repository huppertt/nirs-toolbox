function S = quantizestates(q,S)
%QUANTIZESTATES   

%   Author(s): V. Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.


if strcmpi(class(S),'filtstates.dfiir'),
    S.Numerator = single(double(S.Numerator));
    S.Denominator = single(double(S.Denominator));
else
    S = single(double(S));
end

% [EOF]
