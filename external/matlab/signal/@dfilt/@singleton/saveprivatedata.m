function s = saveprivatedata(this)
%SAVEPRIVATEDATA   Save the private data.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

s = base_saveprivatedata(this);

s.privnormGain = get(this, 'privnormGain');

% [EOF]
