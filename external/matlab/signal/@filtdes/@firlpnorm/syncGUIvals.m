function syncGUIvals(h,arrayh)
%SYNCGUIVALS Sync values from frames.
%
%   Inputs:
%       h - handle to this object   
%       arrayh - array of handles to frames


%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call the super's method
lpnorm_syncGUIvals(h,arrayh);

% Get handles to frames
hf = find(arrayh,'Tag','siggui.filterorder');
hopts = find(arrayh,'Tag','siggui.firlpnormoptsframe');

% Store specs in object
set(h,'order',evaluatevars(get(hf,'order')));

set(h,'initNum',evaluatevars(get(hopts,'InitNum'))); 

set(h,'minPhase',get(hopts,'MinPhase'));
