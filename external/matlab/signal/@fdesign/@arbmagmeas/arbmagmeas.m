function this = arbmagmeas(hfilter, varargin)
%ARBMAGMEAS Construct an ARBMAGMEAS object.

%   Copyright 2005-2011 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));
this = fdesign.arbmagmeas;

normFlag = false;
removeFreqPointFlag = false;

% Parse the inputs.
if length(varargin) > 1
  % Expect specified frequencies to be correctly normalized if a sampling
  % frequency has been specified in the design.
  F = varargin{2};  
  if length(F) < 2
    % Freqz only supports 2 or more frequency points, so add a dummy point
    % instead of returning an error. Remove the point after measuring.
    F = [0 F];
    removeFreqPointFlag = true;
  end
  varargin(2) = [];
  parseconstructorinputs(this, hfilter, varargin{:});
else
  % minfo always returns normalized frequencies
  minfo = parseconstructorinputs(this, hfilter, varargin{:});
  F = minfo.Frequencies;
  normFlag = true; 
end
if this.NormalizedFrequency
  Fs = 2;
else
  Fs = this.Fs;
end

if normFlag
   F = F*Fs/2;
end
this.Frequencies = F;

try
    A = zerophase(hfilter,F,Fs);
catch %#ok<CTCH>
    A = abs(freqz(hfilter,F,Fs));
end

if removeFreqPointFlag
  this.Frequencies = this.Frequencies(2);
  A = A(2);
end

this.Amplitudes = A(:).';

% [EOF]
