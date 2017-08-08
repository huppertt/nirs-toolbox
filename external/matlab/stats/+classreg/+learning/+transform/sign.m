function out = sign(in)

%   Copyright 2010 The MathWorks, Inc.


out = zeros(size(in));
out(in<0) = -1;
out(in>0) = +1;
end
