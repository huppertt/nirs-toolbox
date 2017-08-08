function syncGUIvals(h, arrayh)
%SYNCGUIVALS Sync the specifications from the GUI

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.
    
hopts = find(arrayh, '-class', 'siggui.gremezoptsframe');
set(h, 'ErrorBands', evaluatevars(get(hopts, 'ErrorBands')));

gremez_syncGUIvals(h, arrayh);

% [EOF]
