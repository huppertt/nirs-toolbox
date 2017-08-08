function this = loadobj(this, s)
%LOADOBJ   Load this object

%   Author(s): J. Schickler
%   Copyright 1999-2008 The MathWorks, Inc.

if nargin < 2
    s    = this;
    this = feval(s.class);
end

for indx = 1:length(s.AllSpecs)
    hAllSpecs(indx) = feval(s.AllSpecs{indx}.class);
end

% Backwards compatibility
if isfield(s, 'SpecificationType')
    s.Specification = s.SpecificationType;
end

set(this, ...
    'CapturedState', s.CaptureState, ...
    'AllSpecs',      hAllSpecs);

if strcmpi(this.Specification, s.Specification)
    updatecurrentspecs(this);    
else
    this.Specification = s.Specification;
end

% Allow subclasses to set themselves up before trying to set the specs of
% the contained objects.  If they need that information they can get it
% from S.
thisloadobj(this, s);

% Set up all the properties AFTER setting the specificationtype because
% SETCURRENTSPECS calls SYNCSPECS which would overwrite the specifications.
for indx = 1:length(s.AllSpecs)
    loadobj(this.AllSpecs(indx), s.AllSpecs{indx});
end

% [EOF]
