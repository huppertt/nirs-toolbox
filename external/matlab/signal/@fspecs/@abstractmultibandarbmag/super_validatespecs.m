function [F,E,A,nfpts,Fs,normFreqFlag] = super_validatespecs(this)
%SUPER_VALIDATESPECS   

%   Copyright 2005-2011 The MathWorks, Inc.

% Cache sampling frequency
normFreqFlag = this.privActualNormalizedFrequencyState;
if normFreqFlag
  Fs = 1;
else
  Fs = this.privFs;
end

% Concatenate do care and don't care regions to form frequency and
% amplitudes vectors.
F = this.B1Frequencies;
A = this.B1Amplitudes;
nfpts = length(F);
if nfpts~=length(A),
    error(message('signal:fspecs:abstractmultibandarbmag:super_validatespecs:InvalidSpecifications'))
end
if any(diff(F)<=0)
  error(message('signal:fspecs:abstractmultibandarbmag:super_validatespecs:InvalidSpecificationsUniqueF'))
end

E = [F(1) F(end)];
for i=2:this.NBands,
    nextfband = this.(sprintf('%s%d%s','B',i,'Frequencies'));
    nextA = this.(sprintf('%s%d%s','B',i,'Amplitudes'));
    nfpts = length(nextfband);
    if nfpts~=length(nextA)
      error(message('signal:fspecs:abstractmultibandarbmag:super_validatespecs:InvalidSpecifications'))
    end
    if any(diff(nextfband)<=0)
      % Frequency points in a band must be unique
       error(message('signal:fspecs:abstractmultibandarbmag:super_validatespecs:InvalidSpecificationsUniqueF'))
    end
    
    if nextfband(1)<F(end)
      error(message('signal:fspecs:abstractmultibandarbmag:super_validatespecs:InvalidSpecificationsDontCareRegion'))
    elseif nextfband(1)==F(end) 
      if nextA(1)~=A(end)
        error(message('signal:fspecs:abstractmultibandarbmag:super_validatespecs:InvalidSpecificationsJunction'))   
      end
      prevbandlength  = length(this.(sprintf('%s%d%s','B',i-1,'Frequencies')));
      if length(nextfband) == 1 && prevbandlength == 1
        error(message('signal:fspecs:abstractmultibandarbmag:super_validatespecs:InvalidSpecificationsAdjSingleBands'))              
      end
    end
    
    F = [F nextfband]; %#ok<*AGROW>
    A = [A nextA];
    if nextfband(1)==F(end),
        % Treat adjacent bands as a single do-care region (for firpm
        % designs only)
        E(end) = F(end);
    else
        E = [E nextfband(1) nextfband(end)];
    end
end

% Force row vectors
nfpts = length(F);
F = F(:).';
A = A(:).';
E = E(:).';

% [EOF]
