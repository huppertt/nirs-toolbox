function s = dispstr(Hb, varargin)
%DISPSTR Display string of coefficients.
%   DISPSTR(Hb) returns a string that can be used to display the coefficients
%   of discrete-time filter Hb.
%
%   See also DFILT.   
  
%   Author: P. Costa
%   Copyright 1988-2004 The MathWorks, Inc.
  
Hd = dispatch(Hb);
s = dispstr(Hd, varargin{:});

% [EOF]
