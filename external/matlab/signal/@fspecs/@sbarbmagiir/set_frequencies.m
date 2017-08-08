function frequencies = set_frequencies(this, frequencies)
%SET_FREQUENCIES   PreSet function for the 'frequencies' property.

%   Author(s): V. Pellissier
%   Copyright 2005-2011 The MathWorks, Inc.
%     

if isempty(frequencies), 
    return; 
end

if any(frequencies<0),
    error(message('signal:fspecs:sbarbmagiir:set_frequencies:FreqsMustBePositive'))
end

f1 = frequencies(1);
f2 = frequencies(end);
if ~this.NormalizedFrequency,
   f1 = f1*2/this.Fs; 
   f2 = f2*2/this.Fs; 
end
if abs(f1)>eps || abs(f2 - 1) > eps,
    if this.NormalizedFrequency,
        error(message('signal:fspecs:sbarbmagiir:set_frequencies:InvalidNormFrequencies'))
    else
        error(message('signal:fspecs:sbarbmagiir:set_frequencies:InvalidFrequencies'))
    end
end

if any(diff(frequencies)<0),
    error(message('signal:fspecs:sbarbmagiir:set_frequencies:MustBeIncreasing'))
end



% [EOF]
