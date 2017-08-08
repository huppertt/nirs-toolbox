function zi = ziexpand(Hd,x,zi)
%ZIEXPAND Expand initial conditions for multiple channels when necessary
%   ZI = ZIEXPAND(Hd, X, ZI) 
%
%   This function is intended to only be used by SUPER_FILTER to expand initial
%   conditions. 
%
%   This should be a private method.   

%   Author: Thomas A. Bryan, R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.


error(nargchk(3,3,nargin,'struct'));

[m,ndata] = size(x);
ndata = max(ndata,1);

if ~(isempty(zi) | any(size(zi,2) == [ndata,1])),
	error(message('signal:dfilt:abstractfilter:ziexpand:InvalidDimensions'));
end

if size(zi,2) == 1 && ndata ~= 1,
    % Don't enter if ndata=1 because an unnecessary copy of the states is
    % made
	zi = zi(:,ones(1,ndata));
end

