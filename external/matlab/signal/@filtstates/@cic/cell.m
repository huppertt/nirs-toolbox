function varargout = cell(this)
%CELL   Convert the states object to a cell array.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

validate(this);

[nsections nchannels] = size(this.Integrator);

int  = {};
comb = {};

for indx = 1:nchannels
    for jndx = 1:nsections
        int{end+1}  = this.Integrator(jndx, indx).Value;
        comb{end+1} = this.Comb(jndx, indx).Value;
    end
end

if nargout > 1
    varargout = {int, comb};
else
    varargout = {{int, comb}};
end

% [EOF]
