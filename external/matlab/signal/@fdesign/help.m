%HELP Display help for a particular filter design algorithm.
%   HELP(D, METHOD) display help for the design algorithm METHOD for the
%   current settings of the filter designer D. METHOD must be one of the
%   strings returned by <a href="matlab:help fdesign/designmethods">DESIGNMETHODS</a>.
%
%   % EXAMPLE - Get help for lowpass Butterworth filters.
%   d = fdesign.lowpass;
%   designmethods(d)
%   help(d,'butter')
%
%   See also FDESIGN, FDESIGN/DESIGN, FDESIGN/DESIGNMETHODS,
%   FDESIGN/DESIGNOPTS.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.



% [EOF]
