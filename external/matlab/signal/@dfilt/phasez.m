%PHASEZ Phase response of a discrete-time filter (unwrapped).
%   [Phi,W] = PHASEZ(Hd,N) returns vectors Phi and W containing the phase
%   response of the discrete-time filter Hd, and the frequencies (in
%   radians) at which it is evaluated. The phase response is evaluated at N
%   points equally spaced around the upper half of the unit circle. If you
%   don't specify N, it defaults to 8192.
%
%   If Hd is a vector of filter objects, Phi becomes a matrix.  Each column
%   of the matrix corresponds to each filter in the vector.  If a row
%   vector of frequency points is specified, each row of the matrix
%   corresponds to each filter in the vector.
%
%   PHASEZ(Hd) display the phase response in the Filter Visualization Tool
%   (FVTool).
%
%   For additional parameters, see SIGNAL/PHASEZ.

%   Author: V. Pellissier
%   Copyright 1988-2002 The MathWorks, Inc.

% Help for the p-coded PHASEZ method of DFILT classes.
