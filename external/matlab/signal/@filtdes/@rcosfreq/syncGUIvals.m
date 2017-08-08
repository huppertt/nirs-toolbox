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
set(d, 'TransitionMode', get(fr,'TransitionMode'));

% Sync the cutoff frequency
set(d, 'Fc', evaluatevars(get(fr,'Fc')));

% Sync bandwidth or rolloff
set(d, get(d,'TransitionMode'), evaluatevars(get(fr,get(fr,'TransitionMode'))));

