function s = dispstr(Hd, varargin)
%DISPSTR Display string of coefficients.
%   DISPSTR(Hd) returns a string that can be used to display the coefficients
%   of discrete-time filter Hd.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2005 The MathWorks, Inc.

[lat, lad] = dispstr(Hd.filterquantizer, Hd.privlattice(:), Hd.privLadder(:), varargin{:});

s = char({[getString(message('signal:dfilt:dfilt:Lattice')) ':']
          lat
          [getString(message('signal:dfilt:dfilt:Ladder')) ':']
          lad
         });

% [EOF]
