function [FreqEdgesCell,AmpEdgesCell,Fcell,Acell,NBands] = formatspecs(this)
%FORMATSPECS Format the specs
% FreqEdgesCell - cell array of band edges
% AmpEdgesCell  - cell array of amplitude values at band edges
% Acell         - cell array with vectors that contain amplitude values per band
% Fcell         - cell array with vectors that contain frequency points per band

%   Copyright 2011 The MathWorks, Inc.

% Get amplitudes and frequencies in the format needed for the design.

NBands = this.NBands;

for i=1:NBands,
  F = this.(sprintf('%s%d%s','B',i,'Frequencies'));
  A = this.(sprintf('%s%d%s','B',i,'Amplitudes'));
    
  Fcell{i} = F; %#ok<*AGROW>
  Acell{i} = A;
    
  if length(F) > 1
    FreqEdgesCell{i} = [F(1) F(end)];
    AmpEdgesCell{i} = [A(1) A(end)];
  else
    FreqEdgesCell{i} = F;
    AmpEdgesCell{i} = A;
  end
end

% [EOF]
