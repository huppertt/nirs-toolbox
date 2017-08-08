function syncGUIvals(h,d,arrayh)
%SYNCGUIVALS Sync values frame.
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

% Make sure to sync the freq spec type first
set(d, 'freqSpecType', get(fr,'freqSpecType'));

% Now sync the actual spec
propname = determine_dynamicprop(d,...
    get(d, 'freqSpecType'),set(d, 'freqSpecType'));
set(d, propname, evaluatevars(get(fr,propname)));


	






