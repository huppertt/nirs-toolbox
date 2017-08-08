function specs = ft_whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


specObjs = get(h,'specobjs');

specs = [];
for n = 1:length(specObjs),
    specs = [specs, whichspecs(specObjs(n))];
end
	
	