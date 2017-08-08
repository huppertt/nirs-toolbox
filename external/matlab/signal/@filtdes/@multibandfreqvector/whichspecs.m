function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call super's method
specs = fv_whichspecs(h);

% Replace default value
indx = find(strcmpi({specs.name},'FrequencyVector'));

specs(indx) = setfield(specs(indx),'defval',[0, 8640, 9600, 12000, 12720, 17520,18480,24000]);

