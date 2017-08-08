function hm = measure(this, Hd, varargin)
%MEASURE   Measure the filter.

%   Copyright 2008 The MathWorks, Inc.

hm = fdesign.bandpassmeas(Hd, this, varargin{:});

% [EOF]
