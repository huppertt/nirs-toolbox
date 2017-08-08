function isnoteq2one = svnoteq2one(q, refsvq)
%SVNOTEQ2ONE Test if Scale values should be treated as wires

%   Returns only non unity scale values in SVQ and their position in the
%   boolean vector ISNOTEQ2ONE

%   Author(s): V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

iseq2one = (refsvq==1);
if isa(refsvq,'embedded.fi'),
    if (refsvq.Signed && refsvq.FractionLength==refsvq.WordLength-1) || ...
            (~refsvq.Signed && refsvq.FractionLength==refsvq.WordLength),
        % Fractional format
        isrealmax = (refsvq==realmax(refsvq));
        iseq2one = iseq2one | isrealmax;
    end
end
isnoteq2one = ~iseq2one;

if isprop(q,'CoeffAutoScale'),
    % Make sure 1 is representable to get fi to round values between
    % 1-eps and 1+eps to 1. (g446098)
    WL = q.CoeffWordLength;
    issigned = q.Signed;
    if issigned, FL = WL-2; else FL = WL-1; end
    svq = fi(refsvq, issigned, WL, FL);
    % Excludes floating-point and fixed-point unit scale values for
    % autoscaling
    isnoteq2one = ~(iseq2one | (svq==1));
end
