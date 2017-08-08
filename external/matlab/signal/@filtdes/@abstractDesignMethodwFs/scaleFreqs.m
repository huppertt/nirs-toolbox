function scaleFreqs(h)
%SCALEFREQS Scale new frequencies according to current Fs and freqUnits.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Do nothing if no properties
if ~isempty(get(h,'dynamicProps')),
    % Get freqUnits
    freqUnits = get(h,'freqUnits');
    
    % Get all possible freqUnits
    freqUnitsOpts = set(h,'freqUnits');
    
    defaultFs = 48000;
    
    switch freqUnits,
    case freqUnitsOpts{1}, % Normalized (0 to 1),
        % Use Fs = 2 for this purpose
        Fs = 2;
    otherwise,
        if ~isempty(findprop(h,'Fs')),
            Fs = get(h,'Fs');
        else
            Fs = defaultFs;
        end
    end
    
    % Get all freq specs
    pfreq = find(get(h,'dynamicProps'),'Description','freqspec');
    
    for n = 1:length(pfreq),
        oldval = get(h,pfreq(n).Name);
        % Scale value
        newval = oldval*Fs/defaultFs;
        set(h,pfreq(n).Name,newval);
    end
end
