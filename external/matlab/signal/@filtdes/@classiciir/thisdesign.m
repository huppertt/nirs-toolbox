function Hd = thisdesign(d)
%THISDESIGN  Local design method.

%   Author(s): R. Losada
%   Copyright 1988-2009 The MathWorks, Inc.


% Frequencies Have been prenormalized (0 to 1)

Hd = design(get(d, 'ResponseTypeSpecs'), d);
    
% Frequencies will be reset to what they were
