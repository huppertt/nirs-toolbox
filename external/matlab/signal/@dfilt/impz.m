%IMPZ Impulse response of discrete-time filter.
%   [H,T] = IMPZ(Hd) computes the impulse response of the discrete-time
%   filter Hd choosing the number of samples for you, and returns the
%   response in column vector H and a vector of times (or sample intervals)
%   in T (T = [0 1 2...]').
%
%   [H,T] = IMPZ(Hd) returns a matrix H if Hd is a vector.  Each column of
%   the matrix corresponds to each filter in the vector. 
%
%   For additional parameters, see SIGNAL/IMPZ.
%
%   See also DFILT, SIGNAL/IMPZ.

%   Author(s): P. Costa
%   Copyright 1988-2002 The MathWorks, Inc.

% Help for the filter's pcoded IMPZ method.

% [EOF]
