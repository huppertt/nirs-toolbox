function s = fin_whichspecs(h)
%FIN_WHICHSPECS Determine which specs are required for this class.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isempty(findtype('iircombTransitionMode')),
    schema.EnumType('iircombTransitionMode', {'Bandwidth', 'Q'});
end

s = cell2struct({'TransitionMode', 'iircombTransitionMode', 'Bandwidth', ...
        {'PropertyPostSet',@transitionmode_listener}, ''}, specfields(h), 2);
s(2) = cell2struct({'Bandwidth','udouble',1200,[],'freqspec'},specfields(h),2);
s(3) = cell2struct({'Q', 'udouble', 35, [], ''}, specfields(h), 2);

% --------------------------------------------------------------------
function transitionmode_listener(hObj, eventData)

d = eventData.AffectedObject;

off = setdiff({'Bandwidth', 'Q'}, eventData.NewValue);

enabdynprop(d, off{:}, 'Off');
enabdynprop(d, eventData.NewValue, 'On');

% [EOF]
