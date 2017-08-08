function dvalue = double(this)
%DOUBLE   Return the double value for the state.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[rows, cols] = size(this);

dvalue = cell(cols, 1);

for indx = 1:rows
    for jndx = 1:cols
        value = get(this(indx, jndx), 'Value');

        if ~isa(value, 'double')
            value = double(value);
        end

        dvalue{jndx} = [dvalue{jndx}; value];
    end
end

dvalue = [dvalue{:}];

% [EOF]
