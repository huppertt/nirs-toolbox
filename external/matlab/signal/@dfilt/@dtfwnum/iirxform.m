function [Ht, anum, aden] = iirxform(Ho,fun,varargin)
%IIRXFORM IIR Transformations

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

% This should be private

[b, a] = tf(Ho);

[num, den, anum, aden] = feval(fun, b, a, varargin{:});

% Create the transformed and allpass filters
Ht  = feval(str2func(class(Ho)), num, den);

% [EOF]
