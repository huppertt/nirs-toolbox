function Apass = minwordlengthApass(f,md,Astop)
%MINWORDLENGTHAPASS Determine the passband ripples of the minimum wordlength filter

% This should be a private method

%   Copyright 2008 The MathWorks, Inc.

% Get Apass from Astop
Rstop = convertmagunits(Astop,'db','linear','stop');
Apass = convertmagunits(Rstop,'linear','db','pass');
