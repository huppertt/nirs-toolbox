function stimcell = defaulttbstimulus(Hb)
%DEFAULTTBSTIMULUS returns a cell array of stimulus types.
%   DEFAULTTBSTIMULUS returns a cell array of stimulus types
%   based on the filter structure of filter object Hq.
%   Possible values are, 'impulse','step','ramp','chirp', and
%   'noise'
%
%   See also DFILT, GENERATETBSTIMULUS.

%   Copyright 2003-2004 The MathWorks, Inc.

stimcell = {'impulse','step','ramp','chirp','noise'};

% [EOF]

