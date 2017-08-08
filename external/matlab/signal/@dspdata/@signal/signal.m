function this = signal(data, fs)
%SIGNAL   Construct a SIGNAL object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

this = dspdata.signal;

set(this, 'Name', 'Signal');

if nargin
    set(this, 'Data', data);
    if nargin > 1
        normalizefreq(this, false, fs);
    end
end

% [EOF]
