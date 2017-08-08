function h = numDenFilterOrder(N,M)
%NUMDENFILTERORDER  Constructor for the num den filter order object.
%
%   Inputs:
%       N - Numerator order
%       M - Denominator order
%
%   Outputs:
%       h - Handle to the design method object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.numDenFilterOrder;


% Set specified values
if nargin > 1, set(h,'numOrder',N); end
if nargin > 1, set(h,'denOrder',N); end







