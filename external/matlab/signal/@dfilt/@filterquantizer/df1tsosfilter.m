function [y,zf] = df1tsosfilter(q,num,den,sv,issvnoteq2one,x,zi);
% DF1TSOSFILTER Filter for DFILT.DF1TSOS class in double precision mode

%   Author(s): V.Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

x = quantizeinput(q,x);

zinum = zi.Numerator;
ziden = zi.Denominator;

% Initialize
y = x;
zfnum = zinum;
zfden = ziden;

% Form 1/a0
den(:,1) = den(:,1).^-1;

% Vectorize num and den (section by section)
num = num.';
num = num(:);
den = den.';
den = den(:);

% Make sure that sv is not empty.  It doesn't matter what we set it to,
% since we only multiply the sv value based on the issvnoteq2one vector.
% This was causing problems when allocating a STL complex object for
% scalevalues (see: mat2cpp in
% toolbox/filterdesign/quanitzation/filters/include/flutilities.h)
if isempty(sv), sv = 0; end

% Make num and den row vectors
num = num.';
den = den.';
[y,zfnum,zfden] = df1tsosfilter(num,den,sv,issvnoteq2one,x,zinum,ziden);

% Interlace num and den states
zf = zi;
zf.Numerator = zfnum;
zf.Denominator = zfden;
