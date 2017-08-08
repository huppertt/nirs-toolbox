function intvalue = int(this)
%INT   Return the integer value of the state.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[rows, cols] = size(this);


intvalue = cell(cols, 1);

for indx = 1:rows
    for jndx = 1:cols
        value = get(this(indx, jndx), 'Value');

        if isa(value, 'double') || isa(value, 'single')

            % If we have a double put it in the largest integer we have, 32 bits.
            value = int32(value);
        elseif isa(value, 'embedded.fi')

            % Use the FI method to convert to the best integer value.
            value = storedInteger(value);
        end


        intvalue{jndx} = [intvalue{jndx}; value];
    end
end

intvalue = [intvalue{:}];

% [EOF]
