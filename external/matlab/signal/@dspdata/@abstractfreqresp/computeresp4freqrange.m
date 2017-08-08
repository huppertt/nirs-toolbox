function [H,W] = computeresp4freqrange(this,userreq_halfrange,ispsd,isnormalized,centerdc) %#ok<INUSD,INUSL>
%COMPUTERESP4FREQRANGE   Calculate the frequency response for the range
%                        requested.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

% Define a boolean flag representing the state of SpectrumRange property.
obj_hashalfrange = ishalfnyqinterval(this);

% Make sure that Fs, frequency, and NormalizedFrequency property are all
% consistent.
normalizefreq(this,logical(isnormalized));

if ~userreq_halfrange && obj_hashalfrange,     % User requested 'whole' but obj has 'half'.
    wholerange(this);

elseif userreq_halfrange && ~obj_hashalfrange, % User requested 'half' but obj has 'whole'.
    halfrange(this);
end
H = this.Data;
W = this.Frequencies;

% [EOF]
