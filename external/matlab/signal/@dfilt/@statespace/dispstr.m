function s = dispstr(Hd,varargin)
%DISPSTR Display string of coefficients.
%   DISPSTR(Hd) returns a string that can be used to display the coefficients
%   of discrete-time filter Hd.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2005 The MathWorks, Inc.
  
[a,b,c,d] = dispstr(Hd.filterquantizer, Hd.privA, Hd.privB, Hd.privC, Hd.privD, varargin{:});

s = char({'A:'
          a
          'B:'
          b
          'C:'
          c
          'D:'
          d
         });

% [EOF]
