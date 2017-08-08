function Ht = firxform(Ho,fun,varargin)
%XFORM Frequency Transformations.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

for n = 1:length(Ho.Stage),   
   Ht(n) = feval(fun, Ho.Stage(n), varargin{:});
end

Ht = feval(str2func(classname(Ho)),Ht(:));

% [EOF]
