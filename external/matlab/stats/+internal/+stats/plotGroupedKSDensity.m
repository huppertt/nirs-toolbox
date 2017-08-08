function varargout = plotGroupedKSDensity(x,group,varargin)
%PLOTGROUPEDKSDENSITY plots grouped kernel smooth densities
%   plotGroupedKSDensity(X,GROUP) plots kernel smooth densities in each
%   group of the data X. X must be a vecotr. GROUP is a grouping variable
%   that can be accepted by GROUP2IDX. Missing data, NaN elments are
%   excluded for computation. No error checking.
%   
%   [XRANGE, PX] = plotGroupedKSDensity(X,GROUP,VARARGIN) computes the HG
%   axis 'xdata' and 'ydata' in each group and plot the lines. XRANGE
%   stores the xdata (a matrix) and PX stores the corresponding ydata (the
%   probablities).
%
%   Optional name/value pairs:
%     Parameter:       Value:
%      'AxisHandle'    A handle to the axis where the plot will be draw on.
%
%      'Color'         A color matrix that specify the color for each
%                      group. The number of rows must match the total
%                      number of groups, nG. The default is line(nG).
%
%      'LineWidth'     A vector that specify the line width property for
%                      each group. Vector length must be nG. The default is
%                      ones(1,nG).
%
%      'LineStyle'     A cell array of strings that specify the line style
%                      property for each group. It must be a vector of
%                      length nG. The default is all solid lines.
%
%      'Width'         Bandwidth of the kernel-smoothing window. The
%                      default is optimal for estimating normal densities,
%                      but you may want to choose a smaller value to reveal
%                      features such as multiple modes. It must be empty or
%                      a vector of length nG.
%
%   Example:
%   load fisheriris
%   [xrang, px] = internal.stats.plotGroupedKSDensity(meas(:,1),species);
%   internal.stats.plotGroupedKSDensity(meas(:,1),species,...
%            'color','bcr','LineWidth',[2 1 2],'LineStyle',{'-',':','-.'});
%   

%   Copyright 2012 The MathWorks, Inc.

if isempty(group)
    group = ones(1,numel(x));
end
grpID = grp2idx(group);

x = x(:);
% Remove missing data
wasNaN = isnan(x)|isnan(grpID);
x(wasNaN) = [];
if isempty(x)
    return;
end
grpID(wasNaN) = [];

grp   = unique(grpID); % unique integer group labels
nGrp  = numel(grp); % total number of groups

paramNames = {'AxisHandle','Color','LineWidth', 'LineStyle', 'Width'};
defaults   = { gca,         hsv(nGrp), ones(1,nGrp)       ,  [], []};

[h,clr,lw,ls, ww] = internal.stats.parseArgs(paramNames,defaults,varargin{:});

cxmax = max(x) ;
cxmin = min(x) ;

if cxmax == cxmin
    [~, xrange]=ksdensity(x);
    xLim = [min(xrange),max(xrange)];
else
    dx = 0.1*range(x) ;
    xLim = [cxmin-dx, cxmax+dx];
    xrange = xLim(1):0.01*dx: xLim(2);
end
px = zeros(nGrp,size(xrange,2));

for i = 1:nGrp
    xg = x(grpID == grp(i) );
    if isempty(ww)
        px(i,:) = ksdensity(xg,xrange);
    else
        px(i,:) = ksdensity(xg,xrange,'Width',ww(i));
    end
end

if isempty(clr)
    clr = lines(nGrp);
elseif ischar(clr) && isvector(clr)
    clr = clr(:);
end

% Now draw the kernel density line of each group
hXLines = plot(h, xrange,px);
% Set the line properties accordingly
for i = 1:nGrp
    set(hXLines(i),'Color',clr(i,:),'LineWidth',lw(i));
    if ~isempty(ls)
        set(hXLines(i),'LineStyle',ls{i});
    end
end
% Draw a black horizontal line as the X axis because axis may be 'off'.
line(xLim,[0 0],'Color','k');

if nargout > 0
    varargout{1} = xrange;
    varargout{2} = px;
end

end
