function staticresponse(this, hax, magunits)
%STATICRESPONSE   

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if this.NormalizedFrequency,
    frequnits = 'normalized (0 to 1)';
else
    frequnits = 'Hz';
end

staticrespengine('setupaxis', hax, frequnits, magunits);

% Allow subclasses to add annotations.
thisstaticresponse(this, hax);

% [EOF]
