function c = dfiltvariable(Hd)
%DFILTVARIABLE DFILT Variable string.

% This should be a private method.

%   Author(s): P. Costa
%   Copyright 1988-2004 The MathWorks, Inc.

c = cell(1,length(Hd.Stage));
for k=1:length(Hd.Stage)
  c{k} = dfiltvariable(Hd.Stage(k));
end

% [EOF]
