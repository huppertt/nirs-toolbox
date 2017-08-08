function [custom, trimmed, warnstr, warnid] = trimcustom(this, Hd)
%TRIMCUSTOM   Trim the custom setting for a given filter.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

custom = get(this, 'UserDefinedSections');
if ~iscell(custom)
    custom = {custom};
end

warnid  = '';
warnstr = '';

% Loop over the custom cell array.  Make sure none of the specified
% indices exceeds NSECTIONS
indx2rm = [];
trimmed = false;
for indx = 1:length(custom)
    exceed = custom{indx} > nsections(Hd);
    if any(exceed)
        warnstr = 'User Defined SOS View exceeds the number of sections.  Ignoring higher section numbers.';
        warnid  = 'exceedsnsecs';
        trimmed = true;
        custom{indx}(exceed) = [];
        if isempty(custom{indx})
            indx2rm = [indx2rm indx];
        end
    end
end
custom(indx2rm) = [];

% [EOF]
