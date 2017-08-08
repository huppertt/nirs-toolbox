function f = base_num(Hb, fcn, varargin)
%BASE_NUM   Gateway for vector support of methods that return numbers.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

f = [];
for indx = 1:length(Hb)
    f = [f feval(fcn, Hb(indx), varargin{:})];
end

% [EOF]
