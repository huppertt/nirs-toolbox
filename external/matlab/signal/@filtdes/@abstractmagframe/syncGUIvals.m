function syncGUIvals(h,d,arrayh)
%SYNCGUIVALS Properties to sync.
%
%   Inputs:
%       h - handle to this object
%		d - handle to design method
%       arrayh - array of handles to frames

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frame to sync from
frames = whichframes(h);

% Find the right frame
fr = find(arrayh,'-class',frames.constructor);

% Make sure to sync the magUnits first
s = mf_whichspecs(h);
set(d, s.name, get(fr,s.name));

% Get the right properties to sync, depending on magUnits
specs = syncspecs(h,d);

% Get properties from frame and set them in the design method
for n = 1:length(specs),
	set(d, specs{n}, evaluatevars(get(fr,specs{n})));
end




