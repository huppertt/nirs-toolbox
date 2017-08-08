function [p, w] = phasez(this, varargin)
%PHASEZ Compute the freqz
%   PHASEZ(H, NFFT, UNITCIRCLE) Compute the freqz for NFFT number of points
%   and a ...

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

if nargin > 1 && isstruct(varargin{end})
    opts = varargin{end};
    varargin(end) = [];
else
    opts.showref  = false;
    opts.showpoly = false;
    opts.sosview  = [];
end

% If there is more than 1 filter, we ignore the sosview settings.
if length(this) > 1 || ~isa(this(1).Filter, 'dfilt.abstractsos')
    opts.sosview = [];
end

[nfft, unitcircle] = freqzparse(this, varargin{:});
[wmin, wmax]       = getfreqinputs(this, unitcircle);

p = {};
w = {};

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

    % If unitcircle == 4, then we want to apply the nfft as frequency
    % points (not as the number of points)
    if length(nfft) > 1,
        [p{end+1}, w{end+1}] = phasez(cFilt, nfft, fs*fsScaling);
        if opts.showref && any(isquantized(cFilt))
            pr = phasez(reffilter(cFilt), nfft, fs*fsScaling);
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

% -------------------------------------------------------------------------
function [p, w] = getresponse(Hd, nfft, same, unitcircle, fs, wmin, wmax)

if unitcircle == 3
    [p, w] = phasez(Hd, linspace(wmin, wmax, nfft), fs);
else

    nfft = max(4,round(nfft*fs/(wmax-wmin)));
    % Get the response of the entire filter normalized.
    p = phasez(Hd, nfft, 'whole');

    if same,
        [p, w] = fastreshape(p, fs, unitcircle, nfft);
    else
        % Complete the response for the filter from wmin to wmax based on
        % the sampling frequency of this individual filter.
        [p, w] = completefreqresp(p, fs, wmin, wmax);
    end
end

% [EOF]
