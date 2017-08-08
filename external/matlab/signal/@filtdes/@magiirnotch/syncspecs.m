function props = syncspecs(h,d)
%SYNCSPECS Properties to sync.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

magUnits = get(d,'magUnits');
magUnitsOpts = set(d,'magUnits');

switch magUnits,
case magUnitsOpts{1}, % 'dB'
	props = {'Apass'};
case magUnitsOpts{2}, % 'Linear'
	props = {'Epass'};
end	

% [EOF]
