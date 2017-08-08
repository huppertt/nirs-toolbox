function varargout = zpk(Hb)
%ZPK  Discrete-time filter zero-pole-gain conversion.
%   [Z,P,K] = ZP(Hb) returns the zeros, poles, and gain corresponding to the
%   discrete-time filter Hb in vectors Z, P, and scalar K respectively.
%
%   See also DFILT.   
  
%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

if length(Hb) > 1,
    error(message('signal:dfilt:basefilter:zpk:InvalidDimensions'));
end

Hd = dispatch(Hb);
[z,p,k] = zpk(Hd);

if nargout
    varargout = {z,p,k};
else
    zplaneplot(z,p);

    title(getString(message('signal:dfilt:dfilt:PoleZeroPlot')));

    % Turn on the grid.
    grid on;
end
