function q = numsamples(x)
%NUMSAMPLES_FAST (STRINGNNDATA)

% Copyright 2010 The MathWorks, Inc.

if isempty(x)
  q = 0;
else
  q = size(x{1,1},2);
end
