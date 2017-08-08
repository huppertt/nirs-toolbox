function [H, W] = freqrespest(this, L, opts, optsstruct)
%FREQRESPEST   Calculate the frequency response estimate.

%   Copyright 1988-2011 The MathWorks, Inc.

if nargin < 4
    optsstruct.showref = true;
    optsstruct.sosview = [];
    if nargin < 3
        opts = dspopts.pseudospectrum;
        opts.NFFT = 512; % NFFT default is 'Nextpow2', but we need a numeric value here.
        if nargin < 2
            L = 10;
        end
    end
end

% If there is more than 1 filter, we ignore the sosview settings.
if length(this) > 1 || ~isa(this(1).Filter, 'dfilt.abstractsos')
    optsstruct.sosview = [];
end

fs = getmaxfs(this);
if isempty(fs)
    fs = 2*pi;
end

if strcmpi(opts.SpectrumRange, 'half')
    wmin = 0;
    wmax = fs/2;
elseif opts.CenterDC
    wmin = -fs/2;
    wmax = fs/2;
else
    wmin = 0;
    wmax = fs;
end

% Always get the data from normalized and twosided.  COMPLETEFREQRESP will
% take care of the rest.
opts.SpectrumRange       = 'whole';
opts.NormalizedFrequency = true;
opts.CenterDC            = false;

nfft = opts.NFFT;

for indx = 1:length(this)
    hindx = this(indx).Filter;
    
    fs = get(this(indx), 'Fs');
    if isempty(fs), fs = getmaxfs(this); end
    if isempty(fs), fs = 2*pi;           end
    
    if ~isempty(optsstruct.sosview)
        hindx = getfilters(optsstruct.sosview, hindx);
    end
    % When sosview set to cumulative, the hindx could be
    % a vector for multiple filters. However, they are either all
    % quantized or all not. So it is safe to call "all".
    if all(isquantized(hindx)) && optsstruct.showref
        hindx = [hindx reffilter(hindx)]; %#ok<*AGROW>
        if isprop(hindx(1),'FromSysObjFlag') && hindx(1).FromSysObjFlag && ...
            ~isempty(hindx(1).ContainedSysObj)
          hindx(2).FromSysObjFlag = true;
          hindx(2).ContainedSysObj = clone(hindx(1).ContainedSysObj);
          release(hindx(2).ContainedSysObj)
        end
    end
    
    for jndx = 1:length(hindx)
        opts.NFFT = force2even(max(4, round(nfft*fs/(wmax-wmin))));
        H{indx}(:, jndx) = freqrespest(hindx(jndx), L, opts);
    end
    H{indx} = convert2db(H{indx});
    
    % Complete the response for the filter from wmin to wmax based on
    % the sampling frequency of this individual filter.
    [H{indx}, W{indx}] = completefreqresp(H{indx}, fs, wmin, wmax);
end

% --------------------------------------------------------
function nfft = force2even(nfft)

if rem(nfft, 2)
    nfft = nfft+1;
end

% [EOF]
