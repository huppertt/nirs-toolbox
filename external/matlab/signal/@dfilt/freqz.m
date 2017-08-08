%FREQZ  Discrete-time filter frequency response.
%   [H,W] = FREQZ(Hd,N) returns the N-point complex frequency response
%   vector H and the N-point frequency vector W in radians/sample of the
%   discrete-time filter Hd.  The frequency response is evaluated at N
%   points equally spaced around the upper half of the unit circle. If N
%   isn't specified, it defaults to 8192.
%
%   FREQZ(Hd) with no output argument will launch FVTool in the Magnitude
%   and Phase Response.
%
%   [H,W] = FREQZ(Hd) returns a matrix H if Hd is a vector.  Each column of
%   the matrix corresponds to each filter in the vector.  If a row vector
%   of frequency points is specified, each row of the matrix corresponds to
%   each filter in the vector.
%
%   For additional parameters, see SIGNAL/FREQZ.
%
%   See also DFILT, SIGNAL/FREQZ, FVTOOL.

% Copyright 1988-2004 The MathWorks, Inc.

% Help for the filter's FREQZ method.

% [EOF]
