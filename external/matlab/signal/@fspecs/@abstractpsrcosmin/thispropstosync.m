function p = thispropstosync(this, p) %#ok<INUSL>
%THISPROPSTOSYNC   

%   Copyright 2008 The MathWorks, Inc.

% Exclude Astop
idx = strmatch('Astop', p);
p(idx) = [];

% [EOF]
