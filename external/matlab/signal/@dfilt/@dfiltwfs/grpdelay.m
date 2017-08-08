function [g, w] = grpdelay(this, varargin)
%GRPDELAY Returns the group delay for the filters

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

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

g  = {};
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
        [g{end+1}, w{end+1}] = grpdelay(cFilt, nfft, fs*fsScaling);
        if opts.showref && any(isquantized(cFilt))
            gr = grpdelay(reffilter(cFilt), nfft, fs*fsScaling);
            g{end} = [g{end} gr];
        end

    else
        inputs = {nfft, same, unitcircle, fs*fsScaling, wmin*fsScaling, wmax*fsScaling};
        [g{end+1}, w{end+1}] = getresponse(cFilt, inputs{:});
        if opts.showref && any(isquantized(cFilt))
            gr = getresponse(reffilter(cFilt), inputs{:});
            g{end} = [g{end} gr];
        end
    end
end

% -------------------------------------------------------------------------
function [g, w] = getresponse(Hd, nfft, same, unitcircle, fs, wmin, wmax)

nfft = max(4,round(nfft*fs/(wmax-wmin)));

g = grpdelay(Hd, nfft, 'whole');

if same,
    [g, w] = fastreshape(g, fs, unitcircle, nfft);
else
    [g, w] = completefreqresp(g, fs, wmin, wmax);
end

% [EOF]
