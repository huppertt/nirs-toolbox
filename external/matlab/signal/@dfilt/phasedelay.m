%PHASEDELAY Phase delay of a digital filter.
%   [PHI,W] = PHASEDELAY(Hd,N) returns the N-point phase delay response
%   vector PHI and the N-point frequency vector W in radians/sample of the
%   filter. The phase response is evaluated at N points equally spaced
%   around the upper half of the unit circle. If N isn't specified, it
%   defaults to 8192.
%
%   If Hd is a vector of filter objects, Phi becomes a matrix.  Each column
%   of the matrix corresponds to each filter in the vector.  If a row
%   vector of frequency points is specified, each row of the matrix
%   corresponds to each filter in the vector.
%
%   PHASEDELAY(Hd) displays the phase delay response in the Filter
%   Visualization Tool (FVTool).
%
%   For additional parameters, see SIGNAL/PHASEDELAY.
%
%   See also DFILT, SIGNAL/PHASEDELAY, FVTOOL.

%   Copyright 1988-2004 The MathWorks, Inc.

% Help for the filter's FREQZ method.
