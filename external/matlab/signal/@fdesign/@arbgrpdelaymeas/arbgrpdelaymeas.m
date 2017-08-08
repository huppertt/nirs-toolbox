function this = arbgrpdelaymeas(hfilter, varargin)
%ARBGRPDELAYMEAS Construct an ARBGRPDELAYMEAS object.

%   Copyright 2010-2011 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));
removeFreqPointFlag = false;
this = fdesign.arbgrpdelaymeas;

% Parse the inputs.
if length(varargin) > 1
  F = varargin{2};
  if length(F) < 2
    % grpdelay only supports 2 or more frequency points, so add a dummy
    % point instead of returning an error. Remove the point after
    % measuring.
    F = [0 F];
    removeFreqPointFlag = true;
  end  
  varargin(2) = [];
  minfo = parseconstructorinputs(this, hfilter, varargin{:});
else
  minfo = parseconstructorinputs(this, hfilter, varargin{:});
  F = minfo.Frequencies;
end

if this.NormalizedFrequency
  Fs = 2;
else
  Fs = this.Fs; 
end
    
% Measure the arbitrary magnitude filter.
this.Frequencies = F;
this.TotalGroupDelay = grpdelay(hfilter,this.Frequencies,Fs);
this.TotalGroupDelay = this.TotalGroupDelay(:).';
this.NomGrpDelay = minfo.NomGrpDelay;
% Convert to seconds if normalized frequency is false
if ~this.NormalizedFrequency
  this.TotalGroupDelay = this.TotalGroupDelay/Fs;
  this.NomGrpDelay = this.NomGrpDelay/Fs;
end

if removeFreqPointFlag
  this.TotalGroupDelay = this.TotalGroupDelay(2);
  this.Frequencies = this.Frequencies(2);
end



% [EOF]
