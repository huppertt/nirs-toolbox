function [z,p,k] = zpk(Hd)
%ZPK  Discrete-time filter zero-pole-gain conversion.
%   [Z,P,K] = ZPK(Hd) returns the zeros, poles, and gain corresponding to the
%   discrete-time filter Hd in vectors Z, P, and scalar K respectively.
%
%   Example:
%     [b,a] = butter(3,.4);
%     Hd = dfilt.df2t(b,a);
%     [z,p,k] = zpk(Hd)
%
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2006 The MathWorks, Inc.

if ~isequal(size(Hd),[1,1]),
    error(message('signal:dfilt:singleton:zpk:InvalidDimensions'));
end         

[b,a] = tf(Hd);
[z,p,k] = tf2zpk(b,a);
