function fr = whichframes(h)
%WHICHFRAMES  Return constructors of frames needed for FDATool.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frame needed by the freqspec
fr.constructor = 'fdadesignpanel.multibandfreqmag';
fr.setops      = {};

