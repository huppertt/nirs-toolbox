function Ht = firxform(Ho,fun,varargin)
%FIRXFORM FIR Transformations

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[b, a] = tf(reffilter(Ho));

num = feval(fun, b, varargin{:});

% Create the transformed filter
Ht = copy(Ho);
arith = Ht.Arithmetic; % Cache setting
Ht.Arithmetic = 'double';
Ht.Numerator = num;
Ht.Arithmetic = arith; % Reset arithmetic


% [EOF]
