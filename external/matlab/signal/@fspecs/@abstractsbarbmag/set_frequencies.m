function frequencies = set_frequencies(this, frequencies)
%SET_FREQUENCIES PreSet function for the 'frequencies' property.

%   Copyright 2005-2012 The MathWorks, Inc.

% Don't do the check if coming from a filter designer object. The check
% will be done at run time.
if this.FromFilterDesigner
  return;
end

if isempty(frequencies), 
    return; 
end
f1 = frequencies(1);
f2 = frequencies(end);
if ~this.NormalizedFrequency,
   f1 = f1*2/this.Fs; 
   f2 = f2*2/this.Fs; 
end

if ((abs(f1)>eps && abs(f1+1)>eps) || abs(f2 - 1) > eps) 
    if this.NormalizedFrequency,
        error(message('signal:fspecs:abstractsbarbmag:set_frequencies:InvalidFrequenciesNormalized'));            
    else
        error(message('signal:fspecs:abstractsbarbmag:set_frequencies:InvalidFrequenciesUnNormalized'));            
    end
end

if any(diff(frequencies)<=0)
    error(message('signal:fspecs:abstractsbarbmag:set_frequencies:InvalidFrequenciesMonotinicallyInc'));
end
