function out = symmetricismax(in)

%   Copyright 2010 The MathWorks, Inc.


out = classreg.learning.transform.ismax(in);

%   Copyright 2010 The MathWorks, Inc.

out(out<0.5) = -1;
end
