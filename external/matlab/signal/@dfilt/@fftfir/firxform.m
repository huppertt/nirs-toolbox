function Ht = firxform(Ho,fun,varargin)
%FIRXFORM FIR Transformations

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[b, a] = tf(Ho);

num = feval(fun, b, varargin{:});

Ht = copy(Ho);
Ht.Numerator = num;

% [EOF]
