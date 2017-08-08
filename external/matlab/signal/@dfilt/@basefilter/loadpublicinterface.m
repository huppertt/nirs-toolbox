function loadpublicinterface(this, s)
%LOADPUBLICINTERFACE   Load the public interface.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% Loop over the field names.  We cannot simply set 'this' with 's' because
% 's' may now contain extra information.
f = fieldnames(this);
for indx = 1:length(f)
    if isfield(s, f{indx})
        set(this, f{indx}, s.(f{indx}));
    end
end

% [EOF]
