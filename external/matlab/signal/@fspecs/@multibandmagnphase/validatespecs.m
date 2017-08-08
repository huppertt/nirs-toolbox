function [N,F,E,H,nfpts] = validatespecs(this)
%VALIDATESPECS Validate the specs

%   Copyright 2005-2010 The MathWorks, Inc.

% Get filter order
N = this.FilterOrder;

% Concatenate do care and don't care regions to form frequency and
% amplitudes vectors.
F = this.B1Frequencies;
H = this.B1FreqResponse;
nfpts = length(F);
if nfpts~=length(H),
    error(message('signal:fspecs:multibandmagnphase:validatespecs:InvalidSpecifications'))
end

E = [F(1) F(end)];
for i=2:this.NBands,
    nextfband = this.(sprintf('%s%d%s','B',i,'Frequencies'));
    nextH = this.(sprintf('%s%d%s','B',i,'FreqResponse'));
    nfpts = length(nextfband);
    if nfpts~=length(nextH)
      error(message('signal:fspecs:multibandmagnphase:validatespecs:InvalidSpecifications'))
    end
    if nextfband(1)<F(end)
      error(message('signal:fspecs:multibandmagnphase:validatespecs:InvalidSpecificationsDontCareRegion'))
    elseif nextfband(1)==F(end) && nextH(1)~=H(end)
      error(message('signal:fspecs:multibandmagnphase:validatespecs:InvalidSpecificationsJunction'))
    end
    F = [F nextfband]; %#ok<*AGROW>
    H = [H nextH];
    if nextfband(1)==F(end),
        % Treat adjacent bands as a single do-care region
        E(end) = F(end);
    else
        E = [E nextfband(1) nextfband(end)];
    end
end

% Force row vectors
nfpts = length(F);
F = F(:).';
H = H(:).';
E = E(:).';

% [EOF]
