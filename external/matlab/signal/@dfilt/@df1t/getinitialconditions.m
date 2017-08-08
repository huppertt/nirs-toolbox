function ic = getinitialconditions(Hd)
%GETINITIALCONDITIONS Get the initial conditions.

%   Copyright 2009 The MathWorks, Inc.

s = Hd.States;
ic.Num = double(s.Numerator);
ic.Den = double(s.Denominator);

% [EOF]
