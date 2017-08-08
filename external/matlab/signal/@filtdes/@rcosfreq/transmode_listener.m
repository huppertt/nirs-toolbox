function transmode_listener(h,eventdata)
%TRANSMODE_LISTENER Callback for listener to the TransitionMode property.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Get handle to design method
d = eventdata.AffectedObject;
transmode = eventdata.NewValue;


% Get all possibilities
transmodeOpts = set(d,'TransitionMode');

switch transmode
case transmodeOpts{1}, %'bandwidth'
    bwcase(h,d);
case transmodeOpts{2}, %'rolloff'
    rolloffcase(h,d);
end

%---------------------------------------------------------------------------
function bwcase(h,d)
% Bandwidth case

templatecase(h,d,1,2);

%---------------------------------------------------------------------------
function rolloffcase(h,d)
% Rolloff case

templatecase(h,d,2,1);

%---------------------------------------------------------------------------
function templatecase(h,d,caseflag,oldcaseflag)
% Handle the repeated code for bandwidth and rolloff cases
% caseflag: 1 - bandwidth
%           2 - rolloff
% oldcaseflag: 1 - bandwidth
%              2 - rolloff

% Get all possibilities
transmodeOpts = set(d,'TransitionMode');

% Normalize frequencies temporarily
freqUnits = get(d,'freqUnits');
freqUnitsOpts = set(d,'freqUnits');
set(d,'freqUnits',freqUnitsOpts{1}); % Normalized (0 to 1)

% Get the current rolloff/bandwidth
R = get(d,transmodeOpts{oldcaseflag}); 

% Convert rolloff/bandwidth into bandwidth/rolloff
fc = get(d,'Fc');
switch caseflag,
case 1, % 'bandwidth'
    R = 2*fc*R;
case 2, % 'rolloff'
    R = R/(2*fc);
end

% Disable the rolloff/bandwidth property
enabdynprop(d,transmodeOpts{oldcaseflag},'off');

% Enable the bandwidth/rolloff and set the value
enabdynprop(d,transmodeOpts{caseflag},'on');
set(d,transmodeOpts{caseflag},R);

% Set frequencies back to what they were
set(d,'freqUnits',freqUnits);
