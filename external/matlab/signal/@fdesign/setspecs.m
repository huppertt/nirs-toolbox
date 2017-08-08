%SETSPECS   Set design specifications on a filter designer.
%   SETSPECS(D, S1, S2, etc.) sets the specifications in the filter
%   designer D in the same way they are specificied in the constructor of
%   D. Use GET(D, 'DESCRIPTION') for a description of S1, S2, etc...
%
%   SETSPECS(D,...,Fs) specifies the sampling frequency (in Hz). In this
%   case, all other frequency specifications are also in Hz.
%
%   SETSPECS(D,...,MAGUNITS) specifies the units for any magnitude
%   specification given. MAGUNITS can be one of the following: 'linear',
%   'dB', or 'squared'. If this argument is omitted, 'dB' is assumed. Note
%   that the magnitude specifications are always converted and stored in dB
%   regardless of how they were specified.
%
%   By default, all frequency specifications are assumed to be in
%   normalized frequency units. Moreover, all magnitude specifications are
%   assumed to be in dB.
%
%   Use SET(D, 'SPECIFICATION') to get the list of all available
%   specification types for the filter designer D.
%
%   % Example #1 - Get the list of specification types for lowpass filters.
%   d = fdesign.lowpass;
%   set(d, 'Specification')
%
%   % Example #2 - Specify normalized frequencies.
%   d = fdesign.lowpass('N,Fc');
%   setspecs(d, 20, .4);
%   d
%
%   % Example #3 - Specify a sampling frequency (in Hz).
%   d = fdesign.lowpass('N,Fc');
%   setspecs(d, 20, 4, 20);
%   d
%
%   % Example #4 - Specify the magnitude units.
%   d = fdesign.lowpass;
%   setspecs(d, .4, .5, .1, .05, 'linear');
%   d
%
%   See also FDESIGN, FDESIGN/NORMALIZEFREQ.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

% [EOF]
