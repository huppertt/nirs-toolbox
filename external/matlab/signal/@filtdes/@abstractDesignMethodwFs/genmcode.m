function b = genmcode(d)
%GENMCODE Generate MCode

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Add Fs
b = sigcodegen.mcodebuffer;

b.addcr(gettitlestr(d));
b.cr;

if isnormalized(d),
    b.addcr('% All frequency values are normalized to 1.');
else
        
    b.addcr('%% All frequency values are in %s.', get(d, 'freqUnits'));
    b.addcr('Fs = %s;  %s', getmcode(d, 'Fs'), '% Sampling Frequency');
end

b.cr;

b.add(thisgenmcode(d));

% [EOF]
