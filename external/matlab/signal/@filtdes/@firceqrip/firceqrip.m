function d = firceqrip
%FIRCEQRIP  Constructor for this design method object.
%
%   Outputs:
%       d - Handle to the design method object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

d = filtdes.firceqrip;

% Call super's constructor
singleOrder_construct(d);


% Set the order value to something appropriate
set(d,'order',100);

% Set the tag
set(d,'Tag','FIR constrained equiripple');

