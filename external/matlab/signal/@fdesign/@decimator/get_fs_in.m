function fs_in = get_fs_in(this, fs_in)
%GET_FS_IN   PreGet function for the 'fs_in' property.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if this.NormalizedFrequency
    fs_in = [];
else
    fs_in = this.Fs;
end

% [EOF]
