function disp(this)
%DISP   Display this object.

%   Copyright 2007 The MathWorks, Inc.

s = get(this);

f = fieldnames(s);
f = {f{4}, f{1}, f{5:6}, f{3}, f{2}, f{7:end}};
s = reorderstructure(s, f{:});

if s.NormalizedFrequency
    s = rmfield(s, {'Fs', 'Fs_in', 'Fs_out'});
end

siguddutils('dispstr', s);

% [EOF]
