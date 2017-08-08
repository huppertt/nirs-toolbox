function c = coefficientvariables(Hd)
%COEFFICIENTVARIABLES 

%   This should be a private method.

%   Author(s): P. Costa
%   Copyright 1988-2004 The MathWorks, Inc.

c = cell(1,length(Hd.Stage));
for k=1:length(Hd.Stage)
  c{k} = coefficientvariables(Hd.Stage(k));
end


% [EOF]
