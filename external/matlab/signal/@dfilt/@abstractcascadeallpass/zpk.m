function [z,p,k] = zpk(this)
%ZPK  Discrete-time filter zero-pole-gain conversion.
%   [Z,P,K] = ZPK(Hd) returns the zeros, poles, and gain corresponding to the
%   discrete-time filter Hd in vectors Z, P, and scalar K respectively.
  
%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

if ~isequal(size(this),[1,1]),
    error(message('signal:dfilt:abstractcascadeallpass:zpk:InvalidDimensions'));
end         

hd = dispatch(this);
[z,p,k] = zpk(hd);

% [EOF]
