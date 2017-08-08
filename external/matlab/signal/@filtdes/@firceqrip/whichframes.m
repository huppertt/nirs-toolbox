function fr = whichframes(h)
%WHICHFRAMES  Return constructors of frames needed for FDATool.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Call the base method
fr = super_whichframes(h);

% Override the options frame for this method
fr(end) = getoptsframe(h.responseTypeSpec);



