function specs = thisgetspecs(this)
%THISGETSPECS Get the specs.

%   Copyright 2010 The MathWorks, Inc.

specs.Frequencies = this.Frequencies;
specs.GroupDelay = this.GroupDelay;
% Nominal group delay
specs.NomGrpDelay = this.NomGrpDelay; % this value is computed at design time

% [EOF]
