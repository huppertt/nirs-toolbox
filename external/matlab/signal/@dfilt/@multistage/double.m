function h = double(this)
%DOUBLE   Cast filter to a double-precision arithmetic version.
%   See help in dfilt/double.

%   Author(s): R. Losada
%   Copyright 2003-2004 The MathWorks, Inc.

% Get array of contained dfilts
secs = this.Stage;

for n = 1:length(secs),
    newsecs(n) = double(secs(n));
end

h = feval(str2func(class(this)),newsecs(:));


% [EOF]
