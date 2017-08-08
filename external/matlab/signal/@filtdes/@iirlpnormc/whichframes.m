function fr = whichframes(h)
%WHICHFRAMES  Return constructors of frames needed for FDATool.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call alternate method
fr = super_whichframes(h);

% Add the options frame
fr(end+1).constructor = 'siggui.iirlpnormcoptsframe';
fr(end).setops        = {};