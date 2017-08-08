function s = dispstr(this, offset)
%DISPSTR   Return the display string.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin < 2
    offset = 6;
end

[nsections, nchannels] = size(this.Integrator);

% nchannels = length(this);

% Show the first channels integrator and comb lengths.
IntState  = sprintf('[%dx%d States]', nsections, nchannels);

[nsections, nchannels] = size(this.Comb);
IntStr='Integrator';
CombStr='Comb';
CombState = sprintf('[%dx%d States]', nsections, nchannels);

s = sprintf('%s: %s\n%s %s: %s', ...
    IntStr,IntState, blanks(offset), CombStr,CombState);

% [EOF]
