function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call super's method
specs = mwv_whichspecs(h);

% Replace default value
indx = find(strcmpi({specs.name},'MagnitudeVector'));

specs(indx) = setfield(specs(indx),'defval',[1 1 0 0 1 1 0 0]);

% Replace default value
windx = find(strcmpi({specs.name},'WeightVector'));

specs(windx) = setfield(specs(windx),'defval',[1 1 1 1]);

