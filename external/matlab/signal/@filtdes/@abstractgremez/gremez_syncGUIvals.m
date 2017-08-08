function gremez_syncGUIvals(h, arrayh)
%GREMEZ_SYNCGUIVALS Sync the GREMEZ specific properties.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isspecify(h),
    hopts = find(arrayh, '-class', 'siggui.gremezoptsframe');
    
    if isdynpropenab(h, 'SinglePointBands'),
        set(h, 'SinglePointBands', evaluatevars(get(hopts, 'SinglePointBands')));
    end
    set(h, 'ForcedFreqPoints', evaluatevars(get(hopts, 'ForcedFreqPoints')));
    if isdynpropenab(h, 'IndeterminateFreqPoints'),
        set(h, 'IndeterminateFreqPoints', evaluatevars(get(hopts, 'IndeterminateFreqPoints')));
    end
else
    wf = 'siggui.gremezoptsframe';
    hopts = find(arrayh, '-class', wf);
    set(h, 'initOrder', evaluatevars(get(hopts, 'initOrder')));
end

set(h, 'DensityFactor', evaluatevars(get(hopts, 'DensityFactor')));
set(h, 'Phase', get(hopts, 'Phase'));
if isdynpropenab(h, 'FIRType'),
    set(h, 'FIRType', get(hopts, 'FIRType'));
end
super_syncGUIvals(h, arrayh);

% [EOF]
