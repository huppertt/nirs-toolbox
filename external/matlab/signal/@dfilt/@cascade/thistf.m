function [num, den] = thistf(Hd)
%THISTF  Convert to transfer function.
%   [NUM,DEN] = THISTF(Hd) converts discrete-time filter Hd to numerator and
%   denominator vectors.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

% The algorithm is equivalent to multiplying rational numbers:
%   b1/a1 * b2/a2 = (b1*b2)/(a1*a2)
% where convolution takes the place of *.

if isempty(Hd.Stage)
  num = [];
  den = [];
  return
end

[num, den] = tf(Hd.Stage(1));

for k=2:length(Hd.Stage)
  [b,a] = tf(Hd.Stage(k));
  num = conv(num,b);
  den = conv(den,a);
end

num = removetrailzeros(num);
den = removetrailzeros(den);
