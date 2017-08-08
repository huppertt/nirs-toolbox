function strs = whichframes(h)
%WHICHFRAMES  Return constructors of frames needed for FDATool.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

c = class(h);

strs.constructor = ['fdadesignpanel' c(findstr(c, '.'):end)];
strs.setops      = {};

% [EOF]
