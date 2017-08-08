function freqUnits_listener(h,eventdata)
%FREQUNITS_LISTENER Callback for listener to the freqUnits property.
    
%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Get oldfreqUnits
oldfreqUnits = get(h,'freqUnits');

if nargin > 1,
    freqUnits = eventdata.NewValue;
    % Get all enabled freqspecs 
    if ~isempty(get(h,'dynamicProps')),
        pfreq = find(get(h,'dynamicProps'),'Description','freqspec','AccessFlags.PublicGet','on',...
            'AccessFlags.PublicSet','on');
    else
        pfreq = [];
    end
else
    freqUnits = oldfreqUnits;
    pfreq = [];
end

% Get all possible freqUnits
freqUnitsOpts = set(h,'freqUnits');


% Enable/disable visibility of Fs depending on units and autoconvert values
switch freqUnits,
case freqUnitsOpts{1}, % Normalized (0 to 1)
    % Get Fs value
    Fs = get(h,'Fs');
    
    % Once we have the value, disable Fs property
    enabdynprop(h,'Fs','off');
    
    normalizevals(h,pfreq,Fs);

case freqUnitsOpts{2}, % Hz
    convertvals(h,pfreq,freqUnitsOpts,oldfreqUnits,@toHz);
    
case freqUnitsOpts{3}, % kHz
    convertvals(h,pfreq,freqUnitsOpts,oldfreqUnits,@tokHz);

case freqUnitsOpts{4}, % MHz
    convertvals(h,pfreq,freqUnitsOpts,oldfreqUnits,@toMHz);
    
case freqUnitsOpts{5}, % GHz
    convertvals(h,pfreq,freqUnitsOpts,oldfreqUnits,@toGHz);
end
    
%-------------------------------------------------------------------------------
function normalizevals(h,pfreq,Fs)
%Normalize frequencies

for n = 1:length(pfreq),
    oldval = get(h,pfreq(n).name);
    newval = 2*oldval/Fs;
    set(h,pfreq(n).name,newval);
end

%-------------------------------------------------------------------------------
function convertvals(h,pfreq,freqUnitsOpts,oldfreqUnits,fcnHndl)
% Convert frequency values

% Enable Fs property
enabdynprop(h,'Fs','on');

% Once enabled, get Fs value
Fs = get(h,'Fs');

if strcmpi(oldfreqUnits,freqUnitsOpts{1}), % Normalized (0 to 1)
    for n = 1:length(pfreq),
        oldval = get(h,pfreq(n).name);
        newval = oldval*Fs/2;
        set(h,pfreq(n).name,newval);
    end 
else
    feval(fcnHndl,h,pfreq,freqUnitsOpts,oldfreqUnits);
end

%--------------------------------------------------------------------------------
function toHz(h,pfreq,freqUnitsOpts,oldfreqUnits)

switch oldfreqUnits,
case freqUnitsOpts{3}, % 'kHz'
    convertbyfactor(h,pfreq,1e3);
case freqUnitsOpts{4}, % 'MHz'
    convertbyfactor(h,pfreq,1e6);
case freqUnitsOpts{5}, % 'GHz'
    convertbyfactor(h,pfreq,1e9);
end

%--------------------------------------------------------------------------------
function tokHz(h,pfreq,freqUnitsOpts,oldfreqUnits)

switch oldfreqUnits,
case freqUnitsOpts{2}, % 'Hz'
    convertbyfactor(h,pfreq,1e-3);
case freqUnitsOpts{4}, % 'MHz'
    convertbyfactor(h,pfreq,1e3);
case freqUnitsOpts{5}, % 'GHz'
    convertbyfactor(h,pfreq,1e6);
end

%--------------------------------------------------------------------------------
function toMHz(h,pfreq,freqUnitsOpts,oldfreqUnits)

switch oldfreqUnits,
case freqUnitsOpts{2}, % 'Hz'
    convertbyfactor(h,pfreq,1e-6);
case freqUnitsOpts{3}, % 'kHz'
    convertbyfactor(h,pfreq,1e-3);
case freqUnitsOpts{5}, % 'GHz'
    convertbyfactor(h,pfreq,1e3);
end

%--------------------------------------------------------------------------------
function toGHz(h,pfreq,freqUnitsOpts,oldfreqUnits)

switch oldfreqUnits,
case freqUnitsOpts{2}, % 'Hz'
    convertbyfactor(h,pfreq,1e-9);
case freqUnitsOpts{3}, % 'kHz'
    convertbyfactor(h,pfreq,1e-6);
case freqUnitsOpts{4}, % 'MHz'
    convertbyfactor(h,pfreq,1e-3);
end

%--------------------------------------------------------------------------------
function convertbyfactor(h,pfreq,g)
% Convert values by given factor.

oldFs = get(h,'Fs');
set(h,'Fs',oldFs*g);
for n = 1:length(pfreq),
    oldval = get(h,pfreq(n).name);
    newval = oldval*g;
    set(h,pfreq(n).name,newval);
end 
    
