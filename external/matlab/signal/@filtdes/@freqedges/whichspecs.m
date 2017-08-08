function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

specs = fv_whichspecs(h);

% Prop name, data type, default value, listener callback
specs(end+1) = cell2struct({'FrequencyEdges','double_vector',...
        [0 9600 12000 24000],[],'freqspec'},specfields(h),2);




