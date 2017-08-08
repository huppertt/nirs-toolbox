function [q, num, den, sv, issvnoteq2one] = dispatchsecfilter(Hd)
%DISPATCHSECFILTER Dispatch info for secfilter

%   Copyright 2008 The MathWorks, Inc.

q = Hd.filterquantizer;
num = Hd.privNum;
den = Hd.privDen;
sv = Hd.privScaleValues;
issvnoteq2one = checksv(Hd);


% [EOF]
