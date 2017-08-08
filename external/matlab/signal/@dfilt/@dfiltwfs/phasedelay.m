function [p, w] = phasedelay(this, varargin)
%PHASEDELAY Phase delay

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if nargin > 1 && isstruct(varargin{end})
    opts = varargin{end};
    varargin(end) = [];
else
    opts.showref  = false;
    opts.showpoly = false;
    opts.sosview  = [];
    opts.normalizedfreq = 'off';
end

% If there is more than 1 filter, we ignore the sosview settings.
if length(this) > 1 || ~isa(this(1).Filter, 'dfilt.abstractsos')
    opts.sosview = [];
end

[nfft, unitcircle] = freqzparse(this, varargin{:});
[wmin, wmax]       = getfreqinputs(this, unitcircle);

p  = {};
w  = {};

same = isfssame(this);

for indx = 1:length(this)
    
    cFilt = this(indx).Filter;

    % If it is a polyphase filter, break it down to multiple filters.
    if opts.showpoly && ispolyphase(cFilt)
        rcFactors = cFilt.getratechangefactors;
        fsScaling = 1/rcFactors(1);
        cFilt = polyphase(cFilt, 'object');
    else
        fsScaling = 1;
    end
    
    if ~isempty(opts.sosview)
        cFilt = getfilters(opts.sosview, cFilt);
    end

    % Get the current sampling frequency.  If the getmaxfs returns [], then
    % all of the filters are normalized.  When this happens we set Fs = 2*pi
    fs = get(this(indx), 'Fs');
    if isempty(fs), fs = getmaxfs(this); end
    if isempty(fs), fs = 2*pi; fsScaling = 1; end

    if unitcircle == 4,
        [p{end+1}, w{end+1}] = phasedelay(cFilt, nfft, fs*fsScaling);
        if opts.showref && any(isquantized(cFilt))
            pr = phasedelay(reffilter(cFilt), nfft, fs*fsScaling);
            p{end} = [p{end} pr];
        end

    else
        inputs = {nfft, same, unitcircle, fs*fsScaling, wmin*fsScaling, wmax*fsScaling};
        [p{end+1}, w{end+1}] = getresponse(cFilt, inputs{:});
        if opts.showref && any(isquantized(cFilt))
            pr = getresponse(reffilter(cFilt), inputs{:});
            p{end} = [p{end} pr];
        end
    end
end

% To handle the case of switching x-axis between "normalized frequency" and
% "frequency" through context menu.  When the filter is designed with a
% sampling frequency, this sampling frequency is remembered through the
% code and used when passing into the public phasedelay function to do the
% calculation.  Therefore, depending on whether the option of normalized
% frequency is on or off, we need to handle the data ourselves here.

if strcmpi(opts.normalizedfreq, 'on')
    maxfs = getmaxfs(this);
    if ~isempty(maxfs),
        for indx = 1:length(this),
            p{indx} = p{indx}/(2*pi)*maxfs;
        end
    end
end

% -------------------------------------------------------------------------
function [p, w] = getresponse(Hd, nfft, same, unitcircle, fs, wmin, wmax)

nfft = max(4,round(nfft*fs/(wmax-wmin)));

p = phasedelay(Hd, nfft, 'whole', fs);

if same,
    [p, w] = fastreshape(p, fs, unitcircle, nfft);
else
    [p, w] = completefreqresp(p, fs, wmin, wmax);
end

% [EOF]
