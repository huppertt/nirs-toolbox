function h = gremezlpmin
%GREMEZLP Construct a GREMEZLP object

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

h = filtdes.gremezbpmin;

filterType_construct(h);

% [EOF]
