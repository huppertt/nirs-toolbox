function [h, w] = freqz(this, varargin)
%FREQZ Compute the freqz
%   Inputs:
%       this    -   The object
%       NFFT    -   The number of points or a freqvector
%       UC      -   {'Half', 'Whole', 'FFTShift'}

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

h = {};
w = {};

same = isfssame(this);

for indx = 1:length(this)
    cFilt = get(this(indx), 'Filter');
    
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
    if length(nfft) > 1
        [h{end+1}, w{end+1}] = freqz(cFilt, nfft, fs*fsScaling);
        if opts.showref && any(isquantized(cFilt))
            hr = freqz(reffilter(cFilt), nfft, fs*fsScaling);
            if size(nfft, 1) == 1,
                h{end} = [h{end}; hr];
            else
                h{end} = [h{end} hr];
            end
        end
    else
        
        inputs = {nfft, same, unitcircle, fs*fsScaling, wmin*fsScaling, wmax*fsScaling};
        [h{end+1} w{end+1}] = getresponse(cFilt, inputs{:});
        
        if opts.showref && any(isquantized(cFilt))
            hr = getresponse(reffilter(cFilt), inputs{:});
            h{end} = [h{end} hr];
        end
    end
end

% -------------------------------------------------------------------------
function [h, w] = getresponse(Hd, nfft, same, unitcircle, fs, wmin, wmax)

nfft = max(4,round(nfft*fs/(wmax-wmin)));

% Get the response of the entire filter normalized.
h = freqz(Hd, nfft, 'whole');

if same,
    
    [h, w] = fastreshape(h, fs, unitcircle, nfft);
else
    % Complete the response for the filter from wmin to wmax based on
    % the sampling frequency of this individual filter.
    [h, w] = completefreqresp(h, fs, wmin, wmax);
end

% [EOF]
