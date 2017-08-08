function DGDF = dg_dfilt(Stages, Label, expandOrientation)
%DG_DFILT Constructor for this class.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,3,nargin,'struct'));

DGDF = filtgraph.dg_dfilt;

DGDF.stage = Stages;

if nargin > 1
    DGDF.label = Label;
end

if nargin > 2
    DGDF.expandOrientation = expandOrientation;
end

