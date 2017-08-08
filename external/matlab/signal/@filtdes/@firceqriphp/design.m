function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Setup common design args
args = setupdesignparams(h,d);

% Get mag specs
[Dstop, Dpass] = getdesignspecs(h, d);

args{3} = [Dstop,Dpass];

b = firceqrip(args{:},'high');

% Construct object
Hd = dfilt.dffir(b);

% [EOF]
