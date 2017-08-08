function str = coeff2str(NL, coef, stage) %#ok<INUSL>
%COEFF2STR  Convert coefficient to string

%    This should be a private method

%   Copyright 1988-2011 The MathWorks, Inc.

narginchk(3,3);

if stage <= length(coef)
    coef1 = coef(stage);
else
    coef1 = 0;
end

str = mat2str(coef1, 18);