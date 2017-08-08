function centerdc = get_centerdc(this, centerdc) %#ok
%GET_CENTERDC   PreGet function for the 'CenterDC' property.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

if ishalfnyqinterval(this),
    centerdc = false;
else
    centerdc = this.privcenterdc;
end

% [EOF]
