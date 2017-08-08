function Hd = maxflat(this, varargin)
%MAXFLAT   Design a FIR maximally flat filter.

%   Copyright 2008 The MathWorks, Inc.

Hd = design(this, 'maxflat', varargin{:});
