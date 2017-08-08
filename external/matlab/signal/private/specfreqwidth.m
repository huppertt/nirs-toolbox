function width = specfreqwidth(W)
%SPECFREQWIDTH Spectral frequency width
%   Obtain an estimate of the width of an arbitrary frequency vector.  
%
%   This file is for internal use only and may be changed in a future 
%   release.
%   
%   Copyright 2014 The MathWorks, Inc.

% perform column vector conversion before checking 2d matrix

% Cast to enforce Precision rules
W = double(W);

% force column vector
W = W(:);

% Determine the width of the rectangle used to approximate the integral.
width = diff(W);

% There are two cases when spectrum is twosided, CenterDC or not.
% In both cases, the frequency samples does not cover the entire
% 2*pi (or Fs) region due to the periodicity.  Therefore, the
% missing freq range has to be compensated in the integral.  The
% missing freq range can be calculated as the difference between
% 2*pi (or Fs) and the actual frequency vector span.  For example,
% considering 1024 points over 2*pi, then frequency vector will be
% [0 2*pi*(1-1/1024)], i.e., the missing freq range is 2*pi/1024.
%
% When CenterDC is true, if the number of points is even, the
% Nyquist point (Fs/2) is exact, therefore, the missing range is at
% the left side, i.e., the beginning of the vector.  If the number
% of points is odd, then the missing freq range is at both ends.
% However, due to the symmetry of the real signal spectrum, it can
% still be considered as if it is missing at the beginning of the
% vector.  Even when the spectrum is asymmetric, since the
% approximation of the integral is close when NFFT is large,
% putting it in the beginning of the vector is still ok.
%
% When CenterDC is false, the missing range is always at the end of
% the frequency vector since the frequency always starts at 0.

% assume a relatively uniform interval
missingWidth = (W(end) - W(1)) / (numel(W) - 1);

% if CenterDC was not specified, the first frequency point will
% be 0 (DC).
centerDC = ~isequal(W(1),0);
if centerDC
    width = [missingWidth; width];
else
    width = [width; missingWidth];
end

% [EOF]