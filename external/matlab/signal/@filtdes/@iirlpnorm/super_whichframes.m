function fr = super_whichframes(h)
%WHICHFRAMES  Return constructors of frames needed for FDATool.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call the base method
fr = dm_whichframes(h);

% Call the num den order method
fr(end+1) = whichframes(get(h,'numDenFilterOrderObj'));




