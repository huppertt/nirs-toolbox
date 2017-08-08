function phases = get_phases(this, ~)
%GET_PHASES PreGet function for the 'phases' property.

%   Copyright 2005-2010 The MathWorks, Inc.

% Linear Phase
phases = angle(this.FreqResponse);


% [EOF]
