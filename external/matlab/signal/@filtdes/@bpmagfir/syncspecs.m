function props = syncspecs(h,d)
%SYNCSPECS Properties to sync.
%
%   Inputs:
%       h - handle to this object
%		d - handle to container objects
%
%   Outputs:
%       props - cell array, with list of properties to sync from
%               and list of properties to sync to.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

magUnits = get(d,'magUnits');
magUnitsOpts = set(d,'magUnits');

switch magUnits,
case magUnitsOpts{1}, % 'dB'
	props = {'Astop1','Apass','Astop2'};
case magUnitsOpts{2}, % 'Linear'
	props = {'Dstop1','Dpass','Dstop2'};
end	




