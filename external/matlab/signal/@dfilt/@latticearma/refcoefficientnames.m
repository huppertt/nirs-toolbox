function rcnames = refcoefficientnames(this)
%REFCOEFFICIENTNAMES   

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.

super_rcnames = abslatticerefcoefficientnames(this);

rcnames = {super_rcnames{:},'refladder'};

% [EOF]
