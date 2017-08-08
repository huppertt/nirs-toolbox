function [N,F,E,Gd,nfpts] = validatespecs(this)
%VALIDATESPECS   Validate the specs

%   Copyright 2010 The MathWorks, Inc.

% Get filter order
N = this.FilterOrder;

% Concatenate do care and don't care regions to form frequency and
% amplitudes vectors.
F = this.B1Frequencies;
Gd = this.B1GroupDelay;
nfpts = length(F);
if nfpts~=length(Gd),
  error(message('signal:fspecs:multibandarbgrpdelay:validatespecs:InvalidSpecifications'))
end

E = [F(1) F(end)];
for i=2:this.NBands,
    nextfband = this.(sprintf('%s%d%s','B',i,'Frequencies'));
    nextGd = this.(sprintf('%s%d%s','B',i,'GroupDelay'));
    nfpts = length(nextfband);
    if nfpts~=length(nextGd),
      error(message('signal:fspecs:multibandarbgrpdelay:validatespecs:InvalidSpecifications'))
    end
    if nextfband(1)<F(end),
      error(message('signal:fspecs:multibandarbgrpdelay:validatespecs:InvalidSpecificationsDontCareRegion'))      
    elseif nextfband(1)==F(end) && nextGd(1)~=Gd(end),
      error(message('signal:fspecs:multibandarbgrpdelay:validatespecs:InvalidSpecificationsJunction'))      
    end
    F = [F nextfband]; %#ok<*AGROW>
    Gd = [Gd nextGd];
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
Gd = Gd(:).';
E = E(:).';

% [EOF]
