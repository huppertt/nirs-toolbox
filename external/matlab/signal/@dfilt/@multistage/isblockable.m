function b = isblockable(this)
%ISBLOCKABLE True if the object supports the block method

%   Copyright 2009-2010 The MathWorks, Inc.

nsections = length(this.Stage); 
b = zeros(nsections,1);
for k=1:nsections, 
   b(k) = isblockable(this.Stage(k));
end 
b = all(b);

% [EOF]
