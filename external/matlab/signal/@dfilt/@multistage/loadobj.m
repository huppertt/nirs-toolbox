function this = loadobj(s)
%LOADOBJ   Load this object.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

% This is the COPY case.
if ~isstruct(s)
    s = saveobj(s);
    
    % Copy the filters
    for indx = 1:length(s.Stage)
        s.Stage(indx) = copy(s.Stage(indx));
    end
end

% Construct the object.
this = feval(s.class, 1);

this.PersistentMemory    = s.PersistentMemory;
this.NumSamplesProcessed = s.NumSamplesProcessed;

% We need to do this last so that the setting of "PersistentMemory" doesn't
% set all of the contained objects as well.
this.Stage = s.Stage;

if isfield(s, 'version') && s.version.number > 2
    loadmetadata(this, s);
end

% [EOF]
