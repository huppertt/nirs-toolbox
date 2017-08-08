function [result, errstr, errorObj] = ishdlable(Hb)
%ISHDLABLE True if HDL can be generated for the filter object.
%   ISHDLABLE(Hd) determines if HDL code generation is supported for the
%   filter object Hd and returns true or false.
%
%   The determination is based on the filter structure and the 
%   arithmetic property of the filter.
%
%   The optional second return value is a string that specifies why HDL
%   could not be generated for the filter object Hd.
%
%   See also DFILT, GENERATEHDL.

%   Copyright 2003 The MathWorks, Inc.

  result = logical(0);
  
  errorObj = message('signal:dfilt:fftfir:ishdlable:HdlNotSupported',class(Hb));
  errstr = getString(errorObj);

% [EOF]

