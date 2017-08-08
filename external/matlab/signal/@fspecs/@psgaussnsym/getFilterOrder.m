function filterOrder = getFilterOrder(this)
%GETFILTERORDER Get the filterOrder.

%   Copyright 2008 The MathWorks, Inc.

% Return the filter order
filterOrder = this.NumberOfSymbols * this.SamplesPerSymbol;

% [EOF]
