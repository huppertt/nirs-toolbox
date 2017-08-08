function class = set_class(this, class) %#ok<INUSL>
%SET_CLASS PreSet function for the 'Class' property

%   Copyright 2009 The MathWorks, Inc.

if ~(isequal(class,1) || isequal(class,2))
    error(message('signal:fspecs:audioweightingwtclass:set_class:InvalidClassValue'));
    
end

% [EOF]
