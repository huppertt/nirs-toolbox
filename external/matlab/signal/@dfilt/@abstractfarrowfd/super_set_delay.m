function delay = super_set_delay(this, delay)
%SUPER_SET_DELAY   PreSet function for the 'delay' property.

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

this.reffracdelay = delay;

% Quantize the fracdelay
quantizefd(this);

delay = [];

% [EOF]
