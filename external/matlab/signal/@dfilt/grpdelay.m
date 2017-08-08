%GRPDELAY Group delay of a discrete-time filter.
%   [Gd,W] = GRPDELAY(Hd,N) returns length N vectors Gd and W containing
%   the group delay of the discrete-time filter Hd, and the frequencies (in
%   radians) at which it is evaluated. Group delay is -d{angle(w)}/dw.  The
%   frequency response is evaluated at N points equally spaced around the
%   upper half of the unit circle.  For an FIR filter where N is a power of
%   two, the computation is done faster using FFTs.  If you don't specify
%   N, it defaults to 8192.
%
%   GRPDELAY(Hd) with no output argument will launch FVTool in Group Delay.
%
%   [Gd,W] = GRPDELAY(Hd) returns a matrix H if Hd is a vector.  Each
%   column of the matrix corresponds to each filter in the vector.  If a
%   row vector of frequency points is specified, each row of the matrix
%   corresponds to each filter in the vector.
%
%   For additional parameters, see SIGNAL/GRPDELAY.
%
%   See also DFILT, SIGNAL/GRPDELAY, FVTOOL.

% Copyright 1988-2004 The MathWorks, Inc.

% Help for the filter's GRPDELAY method.

% [EOF]
