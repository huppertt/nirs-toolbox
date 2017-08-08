function validaterefcoeffs(~, prop, val)
%VALIDATEREFCOEFFS   

%   Author(s): V. Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

if ~(strcmpi(class(val), 'double') || strcmpi(class(val), 'single') ...
        || strcmpi(class(val), 'embedded.fi') ...
        || strncmpi(class(val), 'int', 3) || strncmpi(class(val), 'uint', 4)),
    error(message('signal:dfilt:singlefilterquantizer:validaterefcoeffs:invalidDataType', prop, 'fi', 'double', 'int*', 'uint*'));
end

if issparse(val),
    error(message('signal:dfilt:singlefilterquantizer:validaterefcoeffs:invalidSparse', prop));
end

if strcmpi(prop, 'sosMatrix'),
    if any(val(:,4)~=1),
        % a0 must be exactly equal to 1 or an error will occur.
        error(message('signal:dfilt:singlefilterquantizer:validaterefcoeffs:invalidsosMatrix'));
    end
end


% [EOF]
