function h = copy(this)
%COPY   Copy the designer.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

h = feval(class(this));

p = propstocopy(this);
for indx = 1:length(p)
    set(h, p{indx}, get(this, p{indx}));
end

thiscopy(h, this);

% Make sure that we use the old specs
h.AllSpecs = copy(this.AllSpecs);

% Empty out the current specifications so that SYNCSPECS does not change
% our copied specs.
h.CurrentSpecs = [];

if strcmpi(this.Specification, h.Specification)
    updatecurrentspecs(h);    
else
    h.Specification = this.Specification;
end

% [EOF]
