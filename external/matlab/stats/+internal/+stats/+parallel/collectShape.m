function shapeParameters = collectShape(shapeParameters, slice)
%COLLECTSHAPE collects dimensioning info for SMARTFOR sliced outputs.
%
%   COLLECTSHAPE is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.

for i=1:length(slice)
    shapeParameters{i} = size(slice{i});
end

