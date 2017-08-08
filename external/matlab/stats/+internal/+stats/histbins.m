function [centers,edges] = histbins(data,cens,freq,binInfo,F,x)
%HISTBINS Compute bin centers for a histogram
%   [CENTERS,EDGES] = HISTBINS(DATA,CENS,FREQ,BININFO,F,X) computes
%   histogram bin centers and edges for the rule specified in BININFO.  For
%   the Freedman-Diaconis rule, HISTBINS uses the empirical distribution
%   function F evaluated at the values X to compute the IQR.  When there is
%   censoring, HISTBINS cannot compute the Scott rule, and F-D is
%   substituted.
%
%   This function is called by SCATTERHIST with just one argument.


%   Copyright 2001-2014 The MathWorks, Inc.

xmin = min(data);
xmax = max(data);
xrange = xmax - xmin;

if nargin<2
    cens = [];
end
if nargin<3
    freq = [];
end

if isempty(freq)
    n = length(data);
else
    n = sum(freq);
end

if nargin>=4
    rule = binInfo.rule;
else
    rule = 2;
end

% Can't compute the variance for the Scott rule when there is censoring,
% use F-D instead.
if (rule == 2) && ~isempty(cens) && any(cens)
    rule = 1; % Freedman-Diaconis
end

switch rule
case 1 % Freedman-Diaconis
    [centers,edges] = iFreedmanDiaconis(F,x,xrange,xmin,xmax,n);

case 2 % Scott
    [centers,edges] = iScott(data,freq,xmin,xmax,n);

case 3 % number of bins given
    [centers,edges] = iFromNumBins(binInfo,xmin,xmax);
    
case 4 % bins centered on integers
    [centers,edges] = iCenteredOnIntegers(xrange,xmin,xmax);
    
case 5 % bin width given
    [centers,edges] = iFixedWidth(binInfo,xrange,xmin,xmax);
end

end

function [centers,edges] = iFreedmanDiaconis(F,x,xrange,xmin,xmax,n)
% Get "quartiles", which may not actually be the 25th and 75th points
% if there is a great deal of censoring, and compute the IQR.
iqr = diff(interp1q([F;1], [x;x(end)], [.25; .75]));

% Guard against too small an IQR.  This may be because most
% observations are censored, or because there are some extreme
% outliers.
if iqr < xrange ./ 10
    iqr = xrange ./ 10;
end

% Compute the bin width proposed by Freedman and Diaconis, and the
% number of bins needed to span the data.  Use approximately that
% many bins, placed at nice locations.
[centers,edges] = internal.stats.binpicker(xmin, xmax, 'freedmandiaconis', n, iqr);

end

function [centers,edges] = iScott(data,freq,xmin,xmax,n)
idx = isinf(data);
if any(idx)
    data = data(~idx);
    xmin = min(data);
    xmax = max(data);
end
if isempty(freq)
    s = sqrt(var(data));
else
    s = sqrt(var(data,freq));
end

% Compute the bin width proposed by Scott, and the number of bins
% needed to span the data.  Use approximately that many bins,
% placed at nice locations.
[centers,edges] = internal.stats.binpicker(xmin, xmax, 'Scott', n, s);
end

function [centers,edges] = iFromNumBins(binInfo,xmin,xmax)
% Do not create more than 1000 bins.
[centers,edges] = internal.stats.binpicker(xmin, xmax, min(binInfo.nbins,1000));
end

function [centers,edges] = iCenteredOnIntegers(xrange,xmin,xmax)
xscale = max(abs([xmin xmax]));
% If there'd be more than 1000 bins, center them on an appropriate
% power of 10 instead.
if xrange > 1000
    step = 10^ceil(log10(xrange/1000));
    xmin = step*round(xmin/step); % make the edges bin width multiples
    xmax = step*round(xmax/step);
    
% If a bin width of 1 is effectively zero relative to the magnitude of
% the endpoints, use a bigger power of 10.
elseif xscale*eps > 1;
    step = 10^ceil(log10(xscale*eps));  
else
    step = 1;
end
centers = floor(xmin):step:ceil(xmax);
edges = (floor(xmin)-.5*step):step:(ceil(xmax)+.5*step);
end

function [centers,edges] = iFixedWidth(binInfo,xrange,xmin,xmax)
% Do not create more than 1000 bins.
binWidth = max(binInfo.width, xrange/1000);
if (binInfo.placementRule == 1) % automatic placement: anchored at zero
    anchor = 0;
else % anchored
    anchor = binInfo.anchor;
end
leftEdge = anchor + binWidth*floor((xmin-anchor) ./ binWidth);
nbins = max(1,ceil((xmax-leftEdge) ./ binWidth));
edges = leftEdge + (0:nbins) .* binWidth; % get exact multiples
centers = edges(2:end) - 0.5 .* binWidth;
end
