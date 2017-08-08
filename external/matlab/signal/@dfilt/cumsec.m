function cumsec(this)
%CUMSEC   Returns a vector of filters for the cumulative sections.
%   H = CUMSEC(Hd) returns a vector of SOS filter objects with the
%   cumulative sections.
%
%   H = CUMSEC(Hd, INDICES) return a vector of SOS filter objects whose
%   indices into the original filter are in INDICES.
%
%   H = CUMSEC(Hd, INDICES, SECONDARY) uses the secondary scaling points to
%   determine where the sections should be split when SECONDARY is true.
%   SECONDARY is false by default.  This option only has an effect on
%   DF2SOS and DF1TSOS filter objects. For these structures, the secondary
%   scaling points refer to the location between the recursive and the
%   nonrecursive part (i.e. the "middle" of the section).
%
%   CUMSEC(Hd,...) with no output arguments plots the magnitude response of
%   the cumulative sections using FVTOOL.
%
%   % Example: Plot the response of the cumulative sections of a 6th order
%   %          sos filter.
%   Hs = fdesign.lowpass('N,F3db',6,.4);
%   Hd = design(Hs,'butter');
%   H = cumsec(Hd);
%   Hfvt = fvtool(H);
%   legend(Hfvt,'First section','First two sections','Overall filter');
%
%   See also DFILT/SCALE, DFILT/SCALECHECK.

%   Author(s): J. Schickler
%   Copyright 2003-2005 The MathWorks, Inc.



% [EOF]
