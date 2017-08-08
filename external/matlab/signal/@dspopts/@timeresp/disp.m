function disp(this)
%DISP   Display this object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

props = {'NormalizedFrequency'};

if ~this.NormalizedFrequency
    props{end+1} = 'Fs';
end

props{end+1} = 'LengthOption';

if strcmpi(this.LengthOption, 'Specified')
    props{end+1} = 'Length';
end

siguddutils('dispstr', this, props);

% [EOF]
