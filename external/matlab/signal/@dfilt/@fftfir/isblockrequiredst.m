function isblockrequiredst(~)
%ISBLOCKREQUIREDST Check if block method requires a DST license

%   Copyright 2011 The MathWorks, Inc.

[b, ~, ~, errObj] = isspblksinstalled;
if ~b
    error(errObj);
end