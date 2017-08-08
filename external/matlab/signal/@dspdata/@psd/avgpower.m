function pwr = avgpower(this,freqrange)
%AVGPOWER    Average power.
%
%   AVGPOWER is not recommended.
%   Use <a href="matlab:help bandpower">bandpower</a> instead.

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

narginchk(1,2);

% Get the spectrum and frequencies.
Pxx = this.Data;
W   = this.Frequencies;

freqrangespecified = false;
if nargin < 2,
    freqrange = [W(1) W(end)];
else
    if ischar(freqrange) || length(freqrange)~=2 || freqrange(1)<W(1) ||...
            freqrange(2)>W(end),
        error(message('signal:dspdata:psd:avgpower:invalidFrequencyRangeVector', 'FREQRANGE'));
    end
    freqrangespecified = true;
end

% Find indices of freq range requested.
idx1 = find(W<=freqrange(1), 1, 'last');
idx2 = find(W>=freqrange(2), 1);

% Determine the width of the rectangle used to approximate the integral.
width = diff(W);
if freqrangespecified,
    lastRectWidth = 0;  % Don't include last point of PSD data.
    width = [width; lastRectWidth];
else
    % Make sure we include Nyquist and the last point before Fs.
    if strcmpi(this.SpectrumType,'onesided'),
        % Include whole bin width.
        lastRectWidth = width(end);        % Assumes uniform data.
        width = [width; lastRectWidth];
    else
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
        if this.NormalizedFrequency,
            Fs = 2*pi;
        else
            Fs = getfs(this);
        end

        missingWidth = Fs - (W(end)-W(1));


        if this.CenterDC
            width = [missingWidth; width];
        else
            width = [width; missingWidth];
        end
    end
end

% Sum the average power over the range of interest.
pwr = width(idx1:idx2)'*Pxx(idx1:idx2,:);

% [EOF]
