function varargout = checksv(Hd)
%CHECKSV  Check number of scale values
%   CHECKSV(Hd) checks the number of scale values in sos filter Hd.  If the
%   scale values are more than the number of sections plus 1, a warning is
%   generated.
%
%   issvnoteq2one = CHECKSV(Hd) also returns issvnoteq2one which specifies
%   which of the first (section number + 1) scale values are 1s.

%   Copyright 2008-2009 The MathWorks, Inc.

nsecs = Hd.nsections;
issvnoteq2one = Hd.issvnoteq2one;

if length(issvnoteq2one) > nsecs + 1
    warning(message('signal:dfilt:abstractsos:filter:ExtraScaleValues'));
    issvnoteq2one = issvnoteq2one(1:nsecs+1);
end

if nargout
    varargout{1} = issvnoteq2one;
end


% [EOF]
