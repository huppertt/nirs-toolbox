function varargout = thissfcnparams(Hd)
%THISSFCNPARAMS Returns the parameters for SDSPFILTER

% Copyright 1988-2010 The MathWorks, Inc.

error(message('signal:dfilt:notSupported', Hd.FilterStructure));

% [EOF]
