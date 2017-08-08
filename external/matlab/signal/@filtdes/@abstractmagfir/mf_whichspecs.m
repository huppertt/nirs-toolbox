function specs = mf_whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Prop name, data type, default value, listener callback
specs = cell2struct({'magUnits','siggui_magspecs_FIRUnits','dB',{'PropertyPreSet',@magUnits_listener},'filtdes.abstractmagfir'},specfields(h),2);


