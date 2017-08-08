function [d, isfull, type] = designmethods(this, varargin)
%DESIGNMETHODS   Return the design methods for this specification object.

%   Copyright 2008-2013 The MathWorks, Inc.

% Parse the inputs

sigOnlyFlag = false;
isfull      = false;
isminphase  = false;
type        = '';

for indx = 1:length(varargin)
    switch lower(varargin{indx})
        case {'iir', 'fir'}
            type = lower(varargin{indx});
        case 'all'
            type = '';
        case 'full'
            isfull = true;
        case {'minphase', 'minimumphase'}
            isminphase = true;
        case 'signalonly'
            % Got design methods that only require SPT
            sigOnlyFlag = true;
    end
end

if sigOnlyFlag
    alldesigns = getdesignobj(this,[],true);
else
    alldesigns = getdesignobj(this);
end

if isempty(alldesigns)
    d = {};
    return;
end

alldesigns = fieldnames(alldesigns);

alldesigns = setdiff(alldesigns, hiddendesigns(this));

% Hard code the list of design methods so that we can use GETDESIGNOBJ all
% over the place.
if isminphase
    allfirs = {'equiripple'};
else
    allfirs = {'equiripple', 'firls', 'fircls', 'freqsamp', 'ifir', 'kaiserwin', ...
        'lagrange', 'multistage', 'window', 'maxflat'};
end
alliirs = {'butter', 'cheby1', 'cheby2', 'ellip', 'iirlpnorm', 'iirls', ...
    'iirlinphase','ansis142','bell41009'};

fir = intersect(allfirs, alldesigns, 'legacy')';
iir = intersect(alliirs, alldesigns, 'legacy')';

% Return the requested method
switch type
    case 'iir'
        d = iir;
    case 'fir'
        d = fir;
    otherwise
        d = [iir; fir];
end

% If the 'full' flag is passed, convert each of the entries to its full
% name.  We have this map in one place by doing this.
if isfull
    for indx = 1:length(d)
        switch lower(d{indx})
            case 'ellip'
                d{indx} = 'Elliptic';
            case 'butter'
                d{indx} = 'Butterworth';
            case 'cheby1'
                d{indx} = 'Chebyshev type I';
            case 'cheby2'
                d{indx} = 'Chebyshev type II';
            case 'kaiserwin'
                d{indx} = 'Kaiser window';
            case {'equiripple', 'window'}
                d{indx} = [upper(d{indx}(1)) d{indx}(2:end)];
            case 'firls'
                d{indx} = 'FIR least-squares';
            case 'iirlinphase'
                d{indx} = 'IIR quasi-linear phase';
            case 'iirls'
                d{indx} = 'IIR least-squares';
            case 'ifir'
                d{indx} = 'Interpolated FIR';
            case 'iirlpnorm'
                d{indx} = 'IIR least p-norm';
            case 'freqsamp'
                d{indx} = 'Frequency sampling';
            case 'lagrange'
                d{indx} = 'Lagrange interpolation';
            case 'multistage'
                d{indx} = 'Multistage equiripple';
            case 'maxflat'
                d{indx} = 'Maximally flat';
            case 'fircls'
                d{indx} = 'FIR constrained least-squares';
            case 'ansis142'
                d{indx} = 'ANSI S1.42 weighting';
            case 'bell41009'
                d{indx} = 'C-Message Bell 41009 weighting';                
        end
    end
end

% [EOF]
