function spfp = convertspecialprops(hObj)
%CONVERTSPECIALPROPS Convert the special frequency properties for design.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Make a cell of 'n's that is 2ce as long as the number of bands.
spfp = repmat({'n'}, 1, 2*nbands(hObj.ResponseTypeSpecs, hObj));

if isdynpropenab(hObj, 'SinglePointBands'),
    spfp = lclformat(spfp, get(hObj, 'SinglePointBands'), 's');
end

spfp = lclformat(spfp, get(hObj, 'ForcedFreqPoints'), 'f');

if isdynpropenab(hObj, 'IndeterminateFreqPoints'),
    spfp = lclformat(spfp, get(hObj, 'IndeterminateFreqPoints'), 'i');
end

% If all the points are 'n' just return [].
if all(strcmpi(spfp, 'n')),
    spfp = [];
end

% -------------------------------------------------------------------------
function spfp = lclformat(spfp, b, c)

for indx = 1:length(b)
    spfp{b(indx)} = c;
end

% [EOF]
