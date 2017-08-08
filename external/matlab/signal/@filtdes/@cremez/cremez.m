function h = cremez
%CREMEZ Construct a CREMEZ object

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.cremez;

% Call super's constructor
singleOrder_construct(h);

% Set the tag
set(h,'Tag','Complex Remez FIR');

% [EOF]
