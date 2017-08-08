function rcnames = refcoefficientnames(this)
%REFCOEFFICIENTNAMES   

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.

super_rcnames = dtfwnumrefcoefficientnames(this);

rcnames = {super_rcnames{:},'refden'};

% [EOF]
