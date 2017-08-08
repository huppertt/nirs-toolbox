function h = iircomb
%IIRCOMB Construct an IIRCOMB object

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

h = filtdes.iircomb;

set(h, 'Tag', 'IIR Comb');

singleOrder_construct(h);

% [EOF]
