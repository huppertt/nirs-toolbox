function varargout = info(this)
%INFO   Returns information about the ps or psd object.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isempty(this.Metadata)
    f = {};
    v = {};
else
    [f, v] = info(this.Metadata);
end

f = {'Type',    f{1}, 'Length',                         getrangepropname(this), 'Fs',    f{2:end}};
v = {this.Name, v{1}, sprintf('%d', length(this.Data)), this.(f{4}),           this.Fs, v{2:end}};

i = cell2struct(v, f, 2);

if nargout
    varargout = {i};
else
    disp(i);
end

% [EOF]
