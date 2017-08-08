function sc = df1df2tscalecheck(Hd,pnorm)
%DF1DF2TSCALECHECK   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

if nargin < 2,
    pnorm = 'Linf';
end

sc = cumnorm(Hd,pnorm,true);


% [EOF]
