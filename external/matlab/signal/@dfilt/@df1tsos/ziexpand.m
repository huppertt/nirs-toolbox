function zi = ziexpand(Hd,x,zi)
%ZIEXPAND Expand initial conditions for multiple channels when necessary
%   ZI = ZIEXPAND(Hd, X, ZI) 
%
%   This function is intended to only be used by SUPER_FILTER to expand initial
%   conditions. 
%
%   This should be a private method.   

%   Author: V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.


error(nargchk(3,3,nargin,'struct'));

[m,ndata] = size(x);
ndata = max(ndata,1);

if ~(isempty(zi) | any(size(zi.Numerator,2) == [ndata,1])),
	error(message('signal:dfilt:df1tsos:ziexpand:InvalidNumeratorStateDimensions'));
end

if ~(isempty(zi) | any(size(zi.Denominator,2) == [ndata,1])),
	error(message('signal:dfilt:df1tsos:ziexpand:InvalidDenominatorStateDimensions'));
end

if size(zi.Numerator,2) == 1,
    zi.Numerator   = zi.Numerator(:,ones(1,ndata));
end
if size(zi.Denominator,2) == 1,
    zi.Denominator = zi.Denominator(:,ones(1,ndata));
end


