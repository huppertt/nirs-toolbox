function [h, w, p] = zerophase(this, varargin)
%ZEROPHASE Compute the zero-phase response
%   ZEROPHASE(H, NFFT, UNITCIRCLE) Compute the zero-phase for NFFT number of points
%   and a ...

%   Author(s): J. Schickler
%   Copyright 1988-2010 The MathWorks, Inc.

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

h  = {};
w  = {};
p  = {};

% Loop over the # of filters
for indx = 1:length(this),
    fs = get(this(indx), 'Fs');
    if isempty(fs), fs = getmaxfs(this); end
    if isempty(fs), fs = 2*pi;           end

    cFilt = this(indx).Filter;

    % If it is a polyphase filter, break it down to multiple filters.
    if opts.showpoly && ispolyphase(cFilt)
        cFilt = polyphase(cFilt, 'object');
    end

%     if opts.showref && any(isquantized(cFilt))
%         cFilt = [cFilt reffilter(cFilt)];
%     end

    if ~isempty(opts.sosview)
        cFilt = getfilters(opts.sosview, cFilt);
    end

    % If the length of nfft > 1, we must have a freqvec
    if length(nfft) > 1,
        
        [h{end+1}, w{end+1}, p{end+1}] = zerophase(cFilt, nfft, fs); %#ok<*AGROW>
        if opts.showref && any(isquantized(cFilt))
            [hr, ~, pr] = zerophase(reffilter(cFilt), nfft, fs);
            
            
             if size(nfft, 1) == 1,
                h{end} = [h{end}; hr];
                p{end} = [p{end}; pr];
                
            else
                h{end} = [h{end} hr];
                p{end} = [p{end} pr];
            end
            
           
        end

        w{end} = nfft;
    else
        
        inputs = {nfft, unitcircle, fs, wmin, wmax};
        [h{end+1}, w{end+1}, p{end+1}] = getresponse(cFilt, inputs{:});
        if opts.showref && any(isquantized(cFilt))
            [hr, ~, pr] = getresponse(reffilter(cFilt), inputs{:});
            h{end} = [h{end} hr];
            p{end} = [p{end} pr];
        end
    end
end

% -------------------------------------------------------------------------
function [h, w, p] = getresponse(Hd, nfft, ~, fs, wmin, wmax)

nfft = max(4, round(nfft*fs/(wmax-wmin)));

% Get the response of the entire filter normalized.
[h, ~, p, opts] = zerophase(Hd, nfft, 'whole');

opts.shift       = 0;
[h, w]           = completefreqresp(h, fs, wmin, wmax, opts);
opts.periodicity = 2;
opts.flip        = 0;
if iscell(p)
    tempp        = p{end}(~isnan(p));
else
    tempp        = p(~isnan(p));
end
opts.shift       = -tempp(end);
p                = completefreqresp(p, fs, wmin, wmax, opts);

% [EOF]
