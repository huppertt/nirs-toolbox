function sv = getsv(Hd,sv)
%GETSV PreGet function for the scale values

%   Copyright 1988-2003 The MathWorks, Inc.

svq = Hd.privScaleValues;
isnoteq2one = Hd.issvnoteq2one;
sv = ones(length(isnoteq2one),1);
% Insert non unity scale values
sv(isnoteq2one) = double(svq);


