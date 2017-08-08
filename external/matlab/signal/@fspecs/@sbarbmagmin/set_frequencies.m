function frequencies = set_frequencies(this, frequencies)
%SET_FREQUENCIES PreSet function for the 'frequencies' property.

%   Copyright 2011 The MathWorks, Inc.

if isempty(frequencies), 
    return; 
end
f1 = frequencies(1);
f2 = frequencies(end);
if ~this.NormalizedFrequency,
   f1 = f1*2/this.Fs; 
   f2 = f2*2/this.Fs; 
end

% Minimum order arbmag designs are only available for real filter designs
if (abs(f1)>eps) || abs(f2 - 1) > eps,
    if this.NormalizedFrequency
      error(message('signal:fspecs:sbarbmagmin:set_frequencies:InvalidFrequenciesNormalized'));      
    else
      error(message('signal:fspecs:sbarbmagmin:set_frequencies:InvalidFrequenciesUnNormalized'));
    end
end

if any(diff(frequencies)<=0)
  error(message('signal:fspecs:sbarbmagmin:set_frequencies:InvalidFrequenciesMonotinicallyInc'));
end

% [EOF]
