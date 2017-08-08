function frequencies = this_setfrequencies(this, frequencies)
%THIS_SETFREQUENCIES PreSet function for the 'frequencies' property.

%   Copyright 2010 The MathWorks, Inc.

if isempty(frequencies), 
    return; 
end
f1 = frequencies(1);
f2 = frequencies(end);
if ~this.NormalizedFrequency,
   f1 = f1*2/this.Fs; 
   f2 = f2*2/this.Fs; 
end
if (abs(f1)>eps) || abs(f2 - 1) > eps,
    if this.NormalizedFrequency
      error(message('signal:fspecs:sbarbgrpdelay:this_setfrequencies:InvalidNormFrequencies'));      
    else
      error(message('signal:fspecs:sbarbgrpdelay:this_setfrequencies:InvalidFrequencies'));
    end
end

if any(diff(frequencies)<0)
  error(message('signal:fspecs:sbarbgrpdelay:this_setfrequencies:InvalidNonMonotonicallyIncreasingFrequencies'));
end

% [EOF]
