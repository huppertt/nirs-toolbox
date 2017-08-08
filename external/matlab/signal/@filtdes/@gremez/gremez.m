function h = gremez
%GREMEZ Constructor for this design method object.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

h = filtdes.gremez;

% Create the dynamic property InitOrder
schema.prop(h, 'InitOrder', 'MATLAB array');

abstractgremez_construct(h, 'gremezOrderMode');

% Create dynamic properties so that they can be turned off depending on the
% ordermode.  InitOrder for minimum.

set(h, 'Tag', 'Generalized REMEZ FIR');

% [EOF]
