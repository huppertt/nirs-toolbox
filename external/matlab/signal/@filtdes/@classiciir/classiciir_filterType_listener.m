function classiciir_filterType_listener(this, eventData)
%CLASSICIIR_FILTERTYPE_LISTENER   Listener to the FilterType property.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

super_filterType_listener(this);

if isspecify(this),
    enabdynprop(this, 'MatchExactly', 'off');
else
    % Disable order
    enabdynprop(this, 'MatchExactly', 'on');
end

% [EOF]
