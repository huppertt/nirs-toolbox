function disp(this)
%DISP   Display this object.

%   Author(s): J. Schickler
%   Copyright 2005-2008 The MathWorks, Inc.

s = get(this);

% Put the Fs_in and Fs_out after Fs.
f = fieldnames(s);
fs_indx = find(strcmpi('Fs', f));
f = {f{1}, f{4:fs_indx}, f{2}, f{3}, f{fs_indx+1:end}};
s = reorderstructure(s, f{:});

if isfield(s, 'SamplesPerSymbol')
    % Put it before specification
    spec_indx = find(strcmpi('Specification', f));
    sps_indx = find(strcmpi('SamplesPerSymbol', f));
    f = {f{1:spec_indx-1}, f{sps_indx}, f{spec_indx:sps_indx-1}, f{sps_indx+1:end}};
    s = reorderstructure(s, f{:});
end

if s.NormalizedFrequency
    s = rmfield(s, {'Fs', 'Fs_in', 'Fs_out'});
end

siguddutils('dispstr', s);

% [EOF]
