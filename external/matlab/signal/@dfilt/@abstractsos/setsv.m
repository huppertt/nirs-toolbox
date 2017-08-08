function sv = setsv(Hd,sv)
%SETSV PreSet function for the ScaleValues property.

%   Copyright 1988-2009 The MathWorks, Inc.

% Make sure that the scale values are stored as a column.
sv = sv(:);
extrasv = ones(nsections(Hd)-length(sv)+1,1);
sv = [sv; extrasv];

% Set the reference and check datatype
set(Hd, 'refScaleValues', sv);

set(Hd, 'NumAddedSV', length(extrasv));

% Quantize the coefficients
quantizecoeffs(Hd);

% Don't duplicate storage
sv = []; 

clearmetadata(Hd);
