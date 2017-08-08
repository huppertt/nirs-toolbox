function syncGUIvals(h,d,arrayh)
%SYNCGUIVALS Sync values frame.
%
%   Inputs:
%       h - handle to this object
%		d - handle to design method
%       arrayh - array of handles to frames

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Get specs to sync to
specs = whichspecs(h);

% Get frame to sync from
frames = whichframes(h);

% Find the right frame
fr = find(arrayh,'-class',frames.constructor);

% Get properties from frame and set them in the design method
for n = 1:length(specs),
	set(d, specs(n).name, evaluatevars(get(fr,specs(n).name)));
end

	






