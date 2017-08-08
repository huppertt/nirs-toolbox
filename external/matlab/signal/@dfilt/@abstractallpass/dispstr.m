function s = dispstr(this, varargin)
%DISPSTR Display string of coefficients.
%   DISPSTR(Hd) returns a string that can be used to display the coefficients
%   of discrete-time filter Hd.
  
  
%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

c = lcldispstr(this,this.AllpassCoefficients(:), varargin{:});

s = char({getString(message('signal:dfilt:dfilt:AllpassCoefficients'))
          c
          });

% [EOF]
