function hm = measure(this, hfilter, varargin)
%MEASURE   Measure the filter.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

hm = fdesign.lowpassmeas(hfilter, this, varargin{:});

% [EOF]
