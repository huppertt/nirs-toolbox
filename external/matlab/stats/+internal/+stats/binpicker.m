function [centers,edges] = binpicker(xmin, xmax, nbins, nobs, extraArg)
%BINPICKER Generate pleasant bin locations for a histogram.
%   CENTERS = BINPICKER(XMIN,XMAX,NBINS) computes centers for histogram
%   bins spanning the range XMIN to XMAX, with extremes of the bins at
%   locations that are a multiple of 1, 2, 3, or 5 times a power of 10.
%
%   CENTERS = BINPICKER(XMIN,XMAX,'FreedmanDiaconis',N,IQR) uses the Freedman-Diaconis
%   rule for bin width to compute the number of bins.  N is the number of
%   data points, and IQR is the sample interquartile range of the data.
%
%   CENTERS = BINPICKER(XMIN,XMAX,'Scott',N,STD) uses Scott's rule for the
%   bin width to compute the number of bins.  N is the number of data
%   points, and STD is the sample standard deviation of the data.  Scott's
%   rule is appropriate for "normal-like" data.
%
%   CENTERS = BINPICKER(XMIN,XMAX,'Sturges',N) uses Sturges' rule for the
%   number of bins.  N is the number of data points.  Sturges' rule tends
%   to give fewer bins than either F-D or Scott.
%
%   For the Freedman-Diaconis, Scott's, or Sturges' rules, BINPICKER
%   automatically generates "nice" bin locations, where the bin width is 1,
%   2, 3, or 5 times a power of 10, and the bin edges fall on multiples of
%   the bin width.  Thus, the actual number of bins will often differ
%   somewhat from the number defined by the requested rule.
%
%   [CENTERS,EDGES] = BINPICKER(...) also returns the bin edges.

%   References:
%      [1] Freedman, D. and P. Diaconis (1981) "On the histogram as a
%          density estimator: L_2 theory", Zeitschrift fur
%          Wahrscheinlichkeitstheorie und verwandte Gebiete, 57:453-476.
%      [2] Scott, D.W. (1979) "On optimal and data-based histograms",
%          Biometrika, 66:605-610.
%      [3] Sturges, H.A. (1926) "The choice of a class interval",
%          J.Am.Stat.Assoc., 21:65-66.

narginchk( 3, 5 );

if nargin < 4
    nobs = [];
end
if nargin < 5
    extraArg = [];
end

if xmax < xmin
    error(message('stats:binpicker:MaxLessThanMin'));
end

if ischar(nbins)
    % Bin width rule specified
    if nobs < 1
        nbins = 1; % give 1 bin for zero-length data
        rule = 0;
    else
        rulenames = {'freedmandiaconis' 'scott' 'sturges' 'sqrt'};
        [~,rule] = internal.stats.getParamVal(nbins,rulenames,'RULE');
    end
else
    % Number of bins specified
    if nbins < 1 || round(nbins) ~= nbins
        error(message('stats:binpicker:NegativeNumBins'));
    end
    rule = 0;
end

xscale = max(abs(xmin),abs(xmax));
xrange = xmax - xmin;

[rawBinWidth, nbins] = iRawBinWidth(extraArg,nobs,nbins,rule,xrange);

% Make sure the bin width is not effectively zero.  Otherwise it will never
% amount to anything, which is what we knew all along.
rawBinWidth = max(rawBinWidth, eps(xscale));
% it may _still_ be zero, if data are all zeroes


% If the data are not constant, place the bins at "nice" locations
if ~iIsConstant(xrange,xscale)
    % Choose the bin width as a "nice" value.
    binWidth = iNiceBinWidth(rawBinWidth);
    
    % Automatic rule specified
    if rule > 0
        [binWidth,rightEdge,leftEdge,nbinsActual] = iAutoBinEdges(xmin,xmax,binWidth);
        
        % Number of bins specified
    else
        [binWidth,rightEdge,leftEdge,nbinsActual] = iBinEdgesFromNumBins(xmin,xmax,nbins,binWidth);
    end
    
else % the data are nearly constant
    % For automatic rules, use a single bin.
    if rule > 0
        nbins = 1;
    end
    
    % There's no way to know what scale the caller has in mind, just create
    % something simple that covers the data.
    [binWidth,rightEdge,leftEdge,nbinsActual] = iBinsForConstant(xmin,xmax,nbins,xscale);
end

edges = linspace( leftEdge, rightEdge, nbinsActual+1 );
centers = edges(2:end) - 0.5 .* binWidth;
end

function [rawBinWidth, nbins] = iRawBinWidth(extraArg,nobs,nbins,rule,xrange)
switch rule
    case 1 % Freedman-Diaconis rule
        % Use the interquartile range to compute the bin width proposed by
        % Freedman and Diaconis, and the number of bins needed to span the
        % data.  Use approximately that many bins, placed at nice
        % locations.
        iqr = extraArg;
        rawBinWidth = 2*iqr ./ nobs.^(1/3);
        
    case 2 % Scott's rule
        % Compute the bin width proposed by Scott, and the number of bins
        % needed to span the data.  Use approximately that many bins,
        % placed at nice locations.
        s = extraArg;
        rawBinWidth = 3.49*s ./ nobs.^(1/3);
        
    case 3 % Sturges' rule for nbins
        nbins = 1 + log2(nobs);
        rawBinWidth = xrange ./ nbins;
        
    case 4 % Sqrt rule for nbins
        nbins = sqrt(nobs);
        rawBinWidth = xrange ./ nbins;
        
    otherwise % number of bins specified
        rawBinWidth = xrange ./ nbins;
end
end

function binWidth = iNiceBinWidth(rawBinWidth)
% iNiceBinWidth   Choose the bin width as a "nice" value.
powOfTen = 10.^floor(log10(rawBinWidth)); % next lower power of 10
relSize = rawBinWidth / powOfTen; % guaranteed in [1, 10)
if  relSize < 1.5
    binWidth = 1*powOfTen;
elseif relSize < 2.5
    binWidth = 2*powOfTen;
elseif relSize < 4
    binWidth = 3*powOfTen;
elseif relSize < 7.5
    binWidth = 5*powOfTen;
else
    binWidth = 10*powOfTen;
end
end

function [binWidth,rightEdge,leftEdge,nbinsActual] = iBinsForConstant(xmin,xmax,nbins,xscale)
% iBinsForConstant   Choose bins for constant data
%
% There's no way to know what scale the caller has in mind, just create
% something simple that covers the data.
if xscale > realmin(class(xscale))
    % Make the bins cover a unit width, or as small an integer width as
    % possible without the individual bin width being zero relative to
    % xscale.  Put the left edge on an integer or half integer below
    % xmin, with the data in the middle 50% of the bin.  Put the left
    % edge similarly above xmax.
    binRange = max(1, ceil(nbins*eps(xscale)));
    leftEdge = floor(2*(xmin-binRange./4))/2;
    rightEdge = ceil(2*(xmax+binRange./4))/2;
else
    leftEdge = -0.5;
    rightEdge = 0.5;
end
binWidth = (rightEdge - leftEdge) ./ nbins;
nbinsActual = nbins;
end

function [binWidth,rightEdge,leftEdge,nbinsActual] = iAutoBinEdges(xmin,xmax,binWidth)
% iAutoBinEdges   Choose bins for "automatic" rules.
%
% Put the bin edges at multiples of the bin width, covering x.  The
% actual number of bins used may not be exactly equal to the requested
% rule. Always use at least two bins.
leftEdge = min(binWidth*floor(xmin ./ binWidth), xmin);
nbinsActual = max(2, ceil((xmax-leftEdge) ./ binWidth));
rightEdge = max(leftEdge + nbinsActual.*binWidth, xmax);
end

function [binWidth,rightEdge,leftEdge,nbinsActual] = iBinEdgesFromNumBins(xmin,xmax,nbins,binWidth)
% iBinEdgesFromNumBins   Choose bins based on number of bins requested
%
% Put the extreme bin edges at multiples of the bin width, covering x.
% Then recompute the bin width to make the actual number of bins used
% exactly equal to the requested number.
leftEdge = min(binWidth*floor(xmin ./ binWidth), xmin);
rightEdge = max(binWidth*ceil(xmax ./ binWidth), xmax);
binWidth = (rightEdge - leftEdge) ./ nbins;
nbinsActual = nbins;
end

function isConstant = iIsConstant(xrange,xscale)
% iIsConstant   True for constant data
isConstant = xrange <= max(sqrt(eps(xscale)), realmin(class(xscale)));
end
