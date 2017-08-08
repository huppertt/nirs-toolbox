function disp(this)
%DISP   Display this object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

props = {'Name', 'Data', 'NormalizedFrequency'};

if ~this.NormalizedFrequency
    props{end+1} = 'Fs';
end

siguddutils('dispstr', this, props);

% [EOF]
