function fr = dm_whichframes(h)
%WHICHFRAMES  Return constructors of frames needed for FDATool.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frames from the filter type
fr = whichframes(get(h,'responseTypeSpecs'));


