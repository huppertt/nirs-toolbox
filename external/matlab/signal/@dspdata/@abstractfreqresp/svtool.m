function svtool(H)
%SVTOOL   Spectral visualization tool for the dspdata objects.

%   Author(s): P. Pacheco
%   Copyright 1988-2006 The MathWorks, Inc.

% First copy the object to de-couple the plot from the command line.
this = copy(H);

% Create a class-specific response object.
hresp = responseobj(this);
rangeopts = getfreqrangeopts(hresp); % rad/sample or Hz
freqrangeidx = getfreqrange(this);

set(hresp,...
    'FrequencyRange',rangeopts{freqrangeidx},...
    'Name',gettitle(this));

plot(hresp);

% [EOF]
