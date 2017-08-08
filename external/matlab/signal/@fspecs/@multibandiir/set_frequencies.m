function frequencies = set_frequencies(this, frequencies)
%SET_FREQUENCIES   PreSet function for the 'frequencies' property.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

if isempty(frequencies), 
    return; 
end

if any(frequencies<0),
    error(message('signal:fspecs:multibandiir:set_frequencies:expectPositiveFrequencies', 'Frequencies'))
end

if any(diff(frequencies)<0),
    error(message('signal:fspecs:multibandiir:set_frequencies:expectIncreasingFrequencies', 'Frequencies'))
end



% [EOF]
