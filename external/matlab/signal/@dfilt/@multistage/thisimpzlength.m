function len = thisimpzlength(Hd, varargin)
%IMPZLENGTH Length of the impulse response for a digital filter.
%   IMPZLENGTH(Hd) returns the length of the impulse response of 
%   the filter defined by Hd.
%  
%   IMPZLENGTH(Hd,TOL) will specify the tolerance for greater or 
%   less accuracy.  By default, TOL = 5e-5.
%  
%   See also IMPZ.

%   Author: J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

% Convert the filter to a transfer function.
[b, a] = tf(Hd);

% Get the length from the function.
len = impzlength(b,a,varargin{:});

% [EOF]
