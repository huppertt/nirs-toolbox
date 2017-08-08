function [F, A] = getmask(this, fcns, ~, specs)
%GETMASK Get the mask.

%   Copyright 2005-2011 The MathWorks, Inc.

if nargin < 4 || isempty(specs)
    specs = getspecs(this);
end

F = [0 specs.Fpass1 specs.Fpass1 specs.Fpass2 specs.Fpass2 1 ...
    NaN 0 specs.Fstop1 specs.Fstop1 specs.Fstop2 specs.Fstop2 1]*fcns.getfs()/2;

astop1 = fcns.formatastop(specs.Astop1);
astop2 = fcns.formatastop(specs.Astop2);

apass = fcns.formatapass(specs.Apass);

A = [astop1(1:2) apass(1) apass(1) astop2(2) astop2(1) NaN ...
    astop1(4) astop1(3) apass(2) apass(2) astop2(3:4)];

% [EOF]
