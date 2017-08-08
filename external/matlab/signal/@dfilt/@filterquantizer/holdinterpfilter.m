function y = holdinterpfilter(this,L,x,ny,nchans)
%HOLDINTERPFILTER   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

% Quantize input
x = quantizeinput(this,x);
y = zeros(ny,nchans);
y = quantizeinput(this,y);

for i=1:nchans,
    xi = x(:,i);
    y(:,i) = reshape(xi(:,ones(L,1)).',ny,1);
end


% [EOF]
