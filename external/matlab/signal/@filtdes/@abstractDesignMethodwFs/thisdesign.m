function Hd = thisdesign(d)
%THISDESIGN  Local design method.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Frequencies Have been prenormalized (0 to 1)

% Call type specific design
h = get(d,'responseTypeSpecs');
Hd = design(h,d);

% Frequencies will be reset to what they were
