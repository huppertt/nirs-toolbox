function Apass = minwordlengthApass(f,md,Astop)
%MINWORDLENGTHAPASS Determine the passband ripples of the minimum wordlength filter

% This should be a private method

%   Copyright 2008 The MathWorks, Inc.

% Use as reference in case original (double-precision) design doesn't meet specs.
if isprop(f,'Apass'),
    Apass = max(md.Apass,f.Apass);
else
    Apass = md.Apass;
end

% [EOF]
