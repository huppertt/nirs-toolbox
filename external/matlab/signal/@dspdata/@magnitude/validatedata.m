function validatedata(this, data)
%VALIDATEDATA   Validate the data

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin < 2
    data = get(this, 'data');
end

% Check that the data is real.
if any(~isreal(data(:)))
    error(message('signal:dspdata:magnitude:validatedata:invalidComplexData', 'DSPDATA.PSD'));
end

% Check that the data is positive.
if any(data(:) < 0)
    error(message('signal:dspdata:magnitude:validatedata:invalidNegativeData', 'DSPDATA.PSD'));
end

% [EOF]
