function [y, t] = timeresp(this, fcn, varargin)
%TIMERESP Calculates the time response

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

w = warning('off', 'FilterDesign:Qfilt:Overflows');

opts.showref  = false;
opts.showpoly = false;
opts.sosview  = [];
N             = [];
while length(varargin)
    if isstruct(varargin{1})
        opts = varargin{1};
        varargin(1) = [];
    elseif isnumeric(varargin{1})
        N = varargin{1};
        varargin(1) = [];
    else
        error(message('signal:dfilt:dfiltwfs:timeresp:InvalidParam'));
    end
end

% If there is more than 1 filter, we ignore the sosview settings.
if length(this) > 1 || ~isa(this(1).Filter, 'dfilt.abstractsos')
    opts.sosview = [];
end

y  = cell(1, length(this));
t  = cell(1, length(this));

% Loop over the quantized filters
for indx = 1:length(this),

    fs = get(this(indx), 'Fs');
    if isempty(fs), fs = getmaxfs(this); end
    if isempty(fs), fs = [];             end

    if isempty(N),
        if isempty(fs),
            lclN = max(getimpzlength(this,opts));
        else
            lclN = ceil(getmaxtime(this, opts)*fs);
        end
    else,
        lclN = N;
    end
    inputs = {lclN, fs};

    cFilt = this(indx).Filter;

    % If it is a polyphase filter, break it down to multiple filters.
    if opts.showpoly && ispolyphase(cFilt)
        cFilt = polyphase(cFilt, 'object');
    end

    if ~isempty(opts.sosview)
        cFilt = getfilters(opts.sosview, cFilt);
    end

    if opts.showref && any(isquantized(cFilt))
        [yq, t{indx}] = fcn(cFilt, lclN, fs);
        yr            = fcn(reffilter(cFilt), lclN, fs);
        y{indx} = [yq yr];
    else
        [y{indx}, t{indx}] = fcn(cFilt, lclN, fs);
    end
end

warning(w);

% --------------------------------------------------------------
function len = getimpzlength(this, opts)
% For the multiple filter case, return the largest length of the
% impulse response so that we can plot that many points

G = get(this, 'Filter');
if ~iscell(G), G = {G}; end

for j = 1:length(G),
    if opts.showpoly && ispolyphase(G{j})
        G{j} = polyphase(G{j}, 'object');
    end
    len(j) = max(impzlength(G{j}));
end

% --------------------------------------------------------------
function time = getmaxtime(this, opts)

len = getimpzlength(this, opts);
fs  = get(this, 'Fs');
if iscell(fs),
    fs  = [fs{:}];
end

time = len./fs;
time = max(time);

% [EOF]
