function [z,p,k] = zpk(Hd)
%ZPK  Discrete-time filter zero-pole-gain conversion.
%   [Z,P,K] = ZPK(Hd) returns the zeros, poles, and gain corresponding to the
%   discrete-time filter Hd in vectors Z, P, and scalar K respectively.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2006 The MathWorks, Inc.
  
z = [];
p = [];
k = 1;
for n=1:length(Hd.Stage)
  [z1,p1,k1] = zpk(Hd.Stage(n));
  z = [z;z1];
  p = [p;p1];
  k = k*k1;
end
