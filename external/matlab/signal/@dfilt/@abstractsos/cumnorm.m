function sc = cumnorm(Hd,pnorm,secondary)
%CUMNORM   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

L = nsections(Hd);
Hdc = cumsec(Hd,1:L,secondary);

sc = zeros(1,L);

for n = 1:L,
    sc(n) = norm(Hdc(n),pnorm);
end

% [EOF]
