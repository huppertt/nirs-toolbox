function super_syncGUIvals(h, arrayh)
%SYNCGUIVALS Sync values from frames.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Call the base method
base_syncGUIvals(h,arrayh);

% Get handle to frame
hf = find(arrayh,'-isa','siggui.abstractfilterorder');

% Store specs in object
mode = get(hf,'mode');
% If orderMode is enabled, set it, otherwise it is in specify mode
if isdynpropenab(h,'orderMode'),
    set(h,'orderMode',mode);
end
modeOpts = set(hf,'mode');
if strcmpi(mode,modeOpts{1}), % Specify
    set(h,'order',evaluatevars(get(hf,'order')));
end


