function hm = measure(this, Hd, varargin)
%MEASURE   Measure the filter.

%   Copyright 2008 The MathWorks, Inc.

hm = fdesign.bandstopmeas(Hd, this, varargin{:});

% [EOF]
