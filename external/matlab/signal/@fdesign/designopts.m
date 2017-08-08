%DESIGNOPTS   Returns the default design parameters.
%   OPTS = DESIGNOPTS(D, METHOD) returns a structure, OPTS, with the default
%   design parameters used by the design method METHOD. METHOD must be one 
%   of the strings returned by <a href="matlab:help fdesign/designmethods">DESIGNMETHODS</a>. 
%
%   Use HELP(D, METHOD) to get a description of the design parameters. 
%
%   % EXAMPLE - Get the design options for minimum order lowpass
%   %           Butterworth filters.
%   d = fdesign.lowpass;
%   designmethods(d)
%   opts = designopts(d, 'butter')
%   help(d,'butter')
%
%   See also FDESIGN, FDESIGN/DESIGN, FDESIGN/DESIGNMETHODS, FDESIGN/HELP.

%   Copyright 2005-2011 The MathWorks, Inc.

% Help file, no code.

% [EOF]
