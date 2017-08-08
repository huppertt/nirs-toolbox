%STEPZ  Discrete-time filter step response.
%   [H,T] = STEPZ(Hd) returns the step response H of the discrete-time
%   filter Hd. The length of column vector H is computed using IMPZLENGTH.
%   The vector of time samples at which H is evaluated is returned in
%   vector T.
%
%   [H,T] = STEPZ(Hd) returns a matrix H if Hd is a vector.  Each column of
%   the matrix corresponds to each filter in the vector. 
%
%   STEPZ(Hd) display the step response in the Filter Visualization Tool
%   (FVTool).
%
%   For additional parameters, see SIGNAL/STEPZ.
%  
%   See also DFILT/IMPZ, DFILT/FREQZ.

%   Author: P. Costa
%   Copyright 1988-2002 The MathWorks, Inc.

% Help for the p-coded STEPZ method of DFILT classes.
