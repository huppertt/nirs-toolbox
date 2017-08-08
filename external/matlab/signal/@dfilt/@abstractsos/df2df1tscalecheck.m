function sc = df2df1tscalecheck(Hd,pnorm)
%DF2DF1TSCALECHECK   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.


if nargin < 2,
    pnorm = 'Linf';
end

L = nsections(Hd);

if all(Hd.ScaleValues(2:end-1) == 1),
    scalevals = false;
else
    scalevals = true;
end

if scalevals,
    sc = zeros(2,L);
else
    sc = zeros(1,L);
end

sc(1,:) = cumnorm(Hd,pnorm,true);

if scalevals,
    sc(2,:) = cumnorm(Hd,pnorm,false);
end

% [EOF]
