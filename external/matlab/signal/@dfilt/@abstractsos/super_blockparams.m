function s = super_blockparams(Hd)
%SUPER_BLOCKPARAMS   

%   Copyright 2006-2011 The MathWorks, Inc.

s = blockparams(Hd.filterquantizer);

refHd = reffilter(Hd);

sv = get(refHd, 'ScaleValues');
nsecs = Hd.nsections;
if length(sv) > nsecs+1
    warning(message('signal:dfilt:basefilter:block:ExtraScaleValues'));
    sv = sv(1:nsecs+1);
end

s.BiQuadCoeffs = mat2str(refHd.sosMatrix,18);
s.ScaleValues  = mat2str(sv,18);
if Hd.OptimizeScaleValues
    s.OptimizeScaleValues = 'on';
else
    s.OptimizeScaleValues = 'off';
end

% [EOF]
