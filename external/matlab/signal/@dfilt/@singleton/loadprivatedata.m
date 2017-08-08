function loadprivatedata(this, s)
%LOADPRIVATEDATA   Load the private data.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

base_loadprivatedata(this, s);

if s.version.number >= 2 && (~isstruct(s) || isfield(s, 'privnormGain')) 
    set(this, 'privnormGain', s.privnormGain);
end

% [EOF]
