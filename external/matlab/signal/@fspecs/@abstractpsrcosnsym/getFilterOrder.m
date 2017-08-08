function filterOrder = getFilterOrder(this)
%GETFILTERORDER Get the filterOrder.

%   Copyright 2008 The MathWorks, Inc.

% Convert filter or der in symbols to filter order in samples and return
filterOrder = this.NumberOfSymbols * this.SamplesPerSymbol;

% [EOF]
