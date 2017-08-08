function s = filt2struct(this)
%FILT2STRUCT   Return a structure representation of the object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% Add the class (dfilt.cascade, dfilt.parallel, mfilt.cascade).
s.class = class(this);

% Add a field for each stage with that stage's structure represenation.
for indx = 1:nstages(this);
    s.(sprintf('Stage%d', indx)) = filt2struct(this.Stage(indx));
end

% [EOF]
