function d = firls
%FIRLS  Constructor for this design method object.
%
%   Outputs:
%       d - Handle to the design method object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

d = filtdes.firls;

% Call super's constructor
singleOrder_construct(d);

% Set the tag
set(d,'Tag','FIR least-squares');





