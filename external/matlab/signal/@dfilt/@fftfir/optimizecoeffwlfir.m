function [Hbest,mrfflag] = optimizecoeffwlfir(this,varargin) %#ok<STOUT,INUSD>
%OPTIMIZECOEFFWLFIR Optimize coefficient wordlength for FIR filters.
%   This should be a private method.

%   Copyright 2009 The MathWorks, Inc.

error(message('signal:dfilt:fftfir:optimizecoeffwlfir:unsupportedFilterStructure'));

% [EOF]
