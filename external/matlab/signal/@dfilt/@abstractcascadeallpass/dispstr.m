function s = dispstr(this, varargin)
%DISPSTR Display string of coefficients.
%   DISPSTR(Hd) returns a string that can be used to display the coefficients
%   of discrete-time filter Hd.
  
  
%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.


coeffs = this.AllpassCoefficients;

fn = fieldnames(coeffs);
for k = 1:length(fn),
    str = lcldispstr(this,getfield(coeffs,fn{k}), varargin{:});
    c{k} = [fn{k}, ' ',str];
end

s = char({getString(message('signal:dfilt:dfilt:AllpassCoefficients'))
          char(c)
          });

% [EOF]
