function p = describe(this)
%DESCRIBE   Describe the object.

%   Copyright 2008 The MathWorks, Inc.

p      = propstoadd(this);
p(1:2) = [];

for indx = 1:length(p),
    
    % Special case FilterOrder and TransitionWidth.
    switch lower(p{indx})
        case 'numpeaksornotches'
            if isequal(lower(this.CombType),'notch')                 
                d = 'Number of Notches';
            else
                d = 'Number of Peaks';
            end
        case 'bw'
            d = 'Bandwidth at GBW';
        case 'gbw'
            d = 'Bandwidth Gain (dB)';
        case 'shelvingfilterorder'
            d = 'Shelving Filter Order';
    end
    p{indx} = d;
end

% Make sure that we get a column.
p = p(:);

% [EOF]
