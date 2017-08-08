function h = fircls
%FIRCLS Construct a FIRCLS design object

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

h = filtdes.fircls;

singleOrder_construct(h);

set(h, 'Tag', 'FIR constrained least-squares');

% [EOF]
