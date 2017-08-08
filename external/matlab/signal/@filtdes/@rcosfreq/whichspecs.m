function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Prop name, data type, default value, listener callback, description
specs(1) = cell2struct({'Fc','udouble',10800,[],'freqspec'},specfields(h),2);

specs(2) = cell2struct({'TransitionMode','freqrcosTransitionMode','rolloff',...
        {'PropertyPreSet',@transmode_listener},'filtdes.rcosfreq'},specfields(h),2);

specs(3) = cell2struct({'bandwidth','udouble',0.5,...
        [],'freqspec'},specfields(h),2);

specs(4) = cell2struct({'rolloff','double0t1',0.5,...
        [],'rolloff'},specfields(h),2);


