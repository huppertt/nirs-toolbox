function fr = whichframes(h)
%WHICHFRAMES  Return constructors of frames needed for FDATool.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

specObjs = get(h,'specobjs');

for n = 1:length(specObjs),
    fr(n) = whichframes(specObjs(n));  
end

