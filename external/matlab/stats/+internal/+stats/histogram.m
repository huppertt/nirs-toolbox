function [hout,Heights,Edges] = histogram(x,varargin)
% HISTOGRAM Plot a histogram.
%    HISTOGRAM(X) plots a histogram of the data in the vector X.  The bin
%    width and the number of bins is selected according to Scott's rule, which
%    is appropriate for data that are "normal-like."
%
%    HISTOGRAM(X,'PARAM1',val1,'PARAM2',val2,...) specifies one or more of
%    the following name/value pairs. You can specify at most one of the
%    parameters NumBins, BinWidth, Edges, and Rule.
%
%        'NumBins'  An integer N specifying the number of bins.
%
%        'BinWidth' A value BW specifying the bin width.
%
%        'Edges'    A vector V specifying the edges of the bins.
%
%        'Rule'     Defines one of the following rules for determining
%                   histogram bins:
%            'Scott'    The bin width is 3.49*STD(X)*LENGTH(X)^(-1/3). This
%                       is the default, and is a good choice for data that
%                       are "normal-like."
%            'FreedmanDiaconis'   The bin width is 2*IQR(X)*LENGTH(X)^(-1/3),
%                       where IQR is the interquartile range of X.  This is
%                       a good choice for data that have much heavier tails
%                       than a normal distribution.
%            'Sturges'  The number of bins is CEIL(1 + LOG2(LENGTH(X))).
%            'Sqrt'     The number of bins is CEIL(SQRT(LENGTH(X))).
%            'Integers' The bins are centered on integers spanning the
%                       range of X. If MAX(X) - MIN(X) > 1000, the bins are
%                       centered on suitable powers of 10 instead.
%
%        'Type'     Defines how the bar heights are chosen:
%            'Counts'   The bar heights are counts, i.e., they sum to
%                       LENGTH(X).  This is the default.
%            'RelFrequency'  The bar heights are relative frequencies, i.e.,
%                       they sum to 1.
%            'Probability' The bar areas are probabilities, i.e. they sum
%                       to 1.  This is a probability density plot. When
%                       EDGES specifies unequal bin widths, 'Probability'
%                       is the most easily interpreted plot. 
%
%        'Style'    Defines one of the following styles:
%            'Ordinary' The bar heights are not cumulative.  This is the
%                       default.
%            'Cumulative' The bar heights are cumulative. This option is
%                       not supported for the 'Probability' type.
%
%    In all cases, the bins span the range of X, and HISTOGRAM adjusts the
%    bin widths or bin counts slightly so that the first and last edges
%    fall on "nice" values. HISTOGRAM plots a maximum of 1000 bins.
%
%    [H,HEIGHTS,EDGES] = HISTOGRAM(...) returns the handle H to the
%    histogram graphics object, a vector HEIGHTS of bar heights, and a
%    vector EDGES of in edges. Depending on the TYPE argument, the values
%    in the vector HEIGHTS are counts, relative frequencies, or
%    probabilities divided by binwidths.
%
%    See also HIST, HISTC, HIST3, BAR.


%   Copyright 2011 The MathWorks, Inc.

%   References:
%      [1] Freedman, D. and P. Diaconis (1981) "On the histogram as a
%          density estimator: L_2 theory", Zeitschrift fur
%          Wahrscheinlichkeitstheorie und verwandte Gebiete, 57:453-476.
%      [2] Scott, D.W. (1979) "On optimal and data-based histograms",
%          Biometrika, 66:605-610.
%      [3] Sturges, H.A. (1926) "The choice of a class interval",
%          J.Am.Stat.Assoc., 21:65-66.
%
% Undocumented type:
%           'Area'           The bar areas are frequencies, i.e. they sum
%                            to LENGTH(X). 'Cumulative' not supported.

x = iEnsureDoubleVector(x);

if isempty(x)
    if nargout > 0
        hout = [];
        Heights = [];
        Edges = [];
    end
    return
end

xmin = min(x); 
xmax = max(x); 
xrange = xmax - xmin;

n = sum(~isnan(x));

okargs =   {'numbins' 'binwidth' 'edges' 'rule'  'type'   'style' ...
            'display' 'outline'};
defaults = {[]        []         ''      'scott' 'counts' 'ordinary' ...
            true      false};
[NumBins,BinWidth,Edges,Rule,Type,Style,Display,Outline,IsSet,PlotArgs] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});

% Process stylistic arguments
[Cumulative,Type,Display,Outline] = ...
    iProcessStyle(Style,Type,Display,Outline,IsSet);

% Get edges for bins
if IsSet.numbins
    Edges = iEdgesFromNumBins(NumBins,xmin,xmax);
elseif IsSet.binwidth
    Edges = iFixedWidthEdges(xrange,xmin,xmax,BinWidth);
elseif IsSet.edges
    Edges = iParseForEdges(Edges);
else
    Edges = iEdgesFromRule(x,Rule,xrange,n,xmin,xmax);
end

BinWidths = diff(Edges);
counts = iComputeCounts(x,Edges);
Heights = iComputeHeights(Type,n,BinWidths,counts,Cumulative);

if Display
    h = iPlotHistogram(Edges,Heights,Outline,PlotArgs{:});
else
    h = [];
end

if nargout > 0
    hout = h; 
end
end

%============= utilities ===============

function y = iqr(x)
% IQR Compute interquartile range.
n = length(x);
F = ((1:n)'-.5) / n;
y = diff(interp1q(F, sort(x(:)), [.25; .75]));
end

%============= Extracted by JC ===============

function h = iPlotHistogram(edges,heights,outline,varargin)
[ax,args] = axescheck(varargin{:});
if isempty(ax)
    ax = gca;
end
if outline
    edges = edges(:)';
    heights = heights(:)';
    x = [edges(1:end-1);edges(1:end-1);edges(2:end)];
    x = [x(:);edges(end)];
    y = [0,heights(1:end-1);heights;heights];
    y = [y(:);0];
    h = plot(x,y,args{:});
else
    h = bar(ax,edges(1:end-1),heights,1,'histc'); % specify all edges but last
    set(h,'FaceColor',[.75,.85,.95],args{:});
end
end

function heights = iComputeHeights(type,n,binWidths,counts,cumulative)
switch lower(type)
case 'counts'
    % sum(heights) is n for 'Ordinary'
    % heights(end) is n for 'Cumulative'
    heights = counts;
case 'relfrequency'
    % sum(heights) is 1 for 'Ordinary'
    % heights(end) is 1 for 'Cumulative'
    heights = counts ./ n;
case 'area'
    % sum(binWidths.*heights) is n for 'Ordinary'
    % binWidths(end).*heights(end) is n for 'Cumulative'
    heights = (counts./binWidths);
case 'probability'
    % sum(binWidths.*heights) is 1 for 'Ordinary'
    % binWidths(end).*heights(end) is 1 for 'Cumulative'
    heights = (counts./binWidths) / sum(counts);
otherwise
    error(message('stats:internal:histogram:BadType'));
end

if cumulative
    heights = cumsum(heights);
end
end

function counts = iComputeCounts(x,edges)
counts = histc(x,edges); 
counts = counts(:)';
counts = cast(counts,superiorfloat(x,edges)); % respect single
counts(end-1) = sum(counts(end-1:end));
counts = counts(1:end-1);
end

function edges = iParseForEdges(edges)
edges = edges(:)';
if ~isnumeric(edges) || ~all(isfinite(edges)) || ~all(diff(edges)>0) ...
                     || numel(edges)<2
    error(message('stats:internal:histogram:BadEdges'));
end
end

function edges = iFreedmanDiaconisEdges(x,xrange,n,xmin,xmax)
if n > 1
    % Guard against too small an IQR.  This may be because there
    % are some extreme outliers.
    iq = max(iqr(x),xrange/10);
    [~,edges] = internal.stats.binpicker(xmin,xmax,'freedmandiaconis',n,iq);
else
    [~,edges] = internal.stats.binpicker(xmin,xmax,1);
end
end

function edges = iScottEdges(x,n,xmin,xmax)
if n > 1
    [~,edges] = internal.stats.binpicker(xmin,xmax,'scott',n,nanstd(x));
else
    [~,edges] = internal.stats.binpicker(xmin,xmax,1);
end
end

function edges = iFixedWidthEdges(xrange,xmin,xmax,binWidth)
% Do not create more than 1000 bins.
if ~isnumeric(binWidth) || ~isscalar(binWidth) || ~isfinite(binWidth) ...
                        || binWidth<=0
    error(message('stats:internal:histogram:BadWidth'));
end

binWidth = max(binWidth, xrange/1000);
leftEdge = xmin - .5*binWidth;
nbins = max(1,ceil((xmax-leftEdge) ./ binWidth));
edges = leftEdge + (0:nbins) .* binWidth; % get exact multiples
end

function edges = iIntegerEdges(xrange,xmin,xmax)
xscale = max(abs([xmin xmax]));
if xrange > 1000
    % If there'd be more than 1000 bins, center them on an appropriate
    % power of 10 instead.
    step = 10^ceil(log10(xrange/1000));
    xmin = step*round(xmin/step); % make the edges bin width multiples
    xmax = step*round(xmax/step);
elseif eps(xscale) > 1;
    % If a bin width of 1 is effectively zero relative to the magnitude of
    % the endpoints, use a bigger power of 10.
    step = 10^ceil(log10(eps(xscale)));
else
    % Otherwise bins are centered on integers.
    step = 1;
end
edges = (floor(xmin)-.5*step):step:(ceil(xmax)+.5*step);
end

function edges = iEdgesFromRule(x,rule,xrange,n,xmin,xmax)
okvals = {'integers' 'sturges' 'sqrt' 'scott' 'freedmandiaconis'};
rule = internal.stats.getParamVal(rule,okvals,'''Rule''');

switch(rule)
case 'integers'
    edges = iIntegerEdges(xrange,xmin,xmax);
case 'sturges'
    [~,edges] = internal.stats.binpicker(xmin,xmax,'sturges',n);
case 'sqrt'
    [~,edges] = internal.stats.binpicker(xmin,xmax,'sqrt',n);
case 'scott'
    edges = iScottEdges(x,n,xmin,xmax);
case 'freedmandiaconis'
    edges = iFreedmanDiaconisEdges(x,xrange,n,xmin,xmax);
end
end

function edges = iEdgesFromNumBins(bins,xmin,xmax)
if ~internal.stats.isScalarInt(bins,1) || ~isfinite(bins)
    error(message('stats:internal:histogram:BadNumBins'));
end

% Do not create more than 1000 bins.
nbins = min(bins,1000);
[~,edges] = internal.stats.binpicker(xmin,xmax,nbins);
end

function x = iEnsureDoubleVector(x)

if ~isvector(x) || ~isreal(x)
    error(message('stats:internal:histogram:BadX'));
elseif ~isfloat(x)
    x = double(x);
end
end

% ------------------------------------------------
function [Cumulative,Type,Display,Outline] = iProcessStyle(Style,Type,Display,Outline,IsSet)

% Process style/type options
Style = internal.stats.getParamVal(Style,{'ordinary' 'cumulative'},'''Style''');
Cumulative = strcmp(Style,'cumulative');

if strcmpi(Type,'area')  % currently not documented
    Type = 'area';
else
    Type = internal.stats.getParamVal(Type,{'counts' 'relfrequency' 'probability'},'''Type''');
end

if Cumulative && any(strcmpi(Type,{'area' 'probability'}))
    error(message('stats:internal:histogram:NoCumulative',Type));
end

Style = internal.stats.getParamVal(Style,{'ordinary' 'cumulative'},'''Style''');
Cumulative = strcmp(Style,'cumulative');

% Process arguments to control plotting
if ~(isscalar(Display) && islogical(Display))
    Display = internal.stats.getParamVal(Display,{'on' 'off'},'''Display''');
    Display = strcmp(Display,'on');
end
if ~(isscalar(Outline) && islogical(Outline))
    Outline = internal.stats.getParamVal(Outline,{'on' 'off'},'''Outline''');
    Outline = strcmp(Outline,'on');
end

% Process arguments for defining bins
if IsSet.numbins+IsSet.binwidth+IsSet.edges+IsSet.rule>1
    error(message('stats:internal:histogram:ConflictingArgs'));
end
end
