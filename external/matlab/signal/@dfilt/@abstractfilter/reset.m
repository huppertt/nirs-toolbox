function reset(Hm)
%RESET Reset the filter.


%   Author: P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.


Hm.NumSamplesProcessed = 0;

% Default the TapIndex to zeros of the appropriate size
ti = Hm.TapIndex;
if isempty(ti), ti = 0; end
Hm.TapIndex = zeros(size(ti));

Hm.States = [];

Hm.nchannels = [];

% Call thisreset so subclasses can do their thing
thisreset(Hm);
