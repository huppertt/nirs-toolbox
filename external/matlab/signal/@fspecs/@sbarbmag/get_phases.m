function phases = get_phases(this, phases)
%GET_PHASES   PreGet function for the 'phases' property.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

% Linear Phase
P = -this.FilterOrder/2*pi;
phases = P*this.Frequencies;


% [EOF]
