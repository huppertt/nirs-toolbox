function b = isspecmet(this, hfdesign)
%ISSPECMET   True if the object is specmet.

%   Copyright 2009 The MathWorks, Inc.

if nargin < 2
    hfdesign = get(this, 'Specification');
end

minfo = measureinfo(hfdesign);

b = true;
n = length(this.Magnitudes);
if any(minfo.A(1:n)<this.Magnitudes)
    b = false;
end
if any(minfo.A(n+2:end)>this.Magnitudes)
    b = false;
end


% [EOF]
