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

% Make sure to sync the transition mode first
set(d, 'DesignType', get(fr,'DesignType'));


