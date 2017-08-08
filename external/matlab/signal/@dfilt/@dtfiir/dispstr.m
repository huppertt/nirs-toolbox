function s = dispstr(Hd, varargin)
%DISPSTR Display string of coefficients.
%   DISPSTR(Hd) returns a string that can be used to display the coefficients
%   of discrete-time filter Hd.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2005 The MathWorks, Inc.

[num, den] = dispstr(Hd.filterquantizer, Hd.privNum(:), Hd.privDen(:), varargin{:});

s = char({[getString(message('signal:dfilt:dfilt:Numerator')) ':']
          num
          [getString(message('signal:dfilt:dfilt:Denominator')) ':']
          den
         });

% [EOF]
