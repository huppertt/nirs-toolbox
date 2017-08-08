function hg = plotGroupedHist(x,group,varargin)  
%PLOTGROUPEDHIST plots grouped histogram
%   plotGroupedHist(X,GROUP) plots histogram in each group of the data X. 
%   X must be a vector. GROUP is a grouping variable that can be accepted 
%   by GROUP2IDX. Missing data, NaN elments are excluded for computation. 
%   No error checking.
%
%   Optional name/value pairs:
%     Parameter:       Value:
%      'AxisHandle'    A handle to the axis where the plot will be draw on.
%      'Color'         A color matrix that specify the color for each
%                      group. The number of rows must match the total
%                      number of groups, nG. The default is line(nG).
%      'AxisOn'        A logic vale to determine whether the axis is 'on'
%                      or not. If 'AxisOn' is false, a black horizontal
%                      line will be drawn as the X axis. The default is
%                      true.
%      'NBins'         A positive integer value, a two-element vector, a
%                      1-by-1, 1-by-2, or a 2-by-1 cell array, specifying
%                      the number of bins for each group in the X and Y
%                      histograms.  All numbers should be positive integers
%                      greater than or equal to 2.
%      'PlotGroup'     A logical value indicating if grouped histograms or
%                      grouped kernel density plots are created when a
%                      grouping variable is given. The default is true. The
%                      histogram of the whole data set will be created if
%                      it is set to false.
%
%   H = PLOTGROUPEDHIST(...) returns an array of handles to the histograms
%   of each group.
%
%   Example:
%   load fisheriris
%   internal.stats.plotGroupedHist(meas(:,1),species,'color','bcr');
%   hg = internal.stats.plotGroupedHist(meas(:,1),species);

%   Copyright 2014 The MathWorks, Inc.

if isempty(group)
    group = ones(1,numel(x));
end
[grpID,gname] = grp2idx(group);

x = x(:);
% Remove missing data
wasNaN = isnan(grpID);
x(wasNaN) = [];
if isempty(x)
    return;
end
grpID(wasNaN) = [];
grp = unique(grpID); % unique integer group labels
nGrp = numel(grp); % total number of groups

paramNames = {'AxisHandle', 'Color', 'AxisOn', 'NBins', 'PlotGroup'};
defaults   = {gca, hsv(nGrp), true, [], true};

[h,clr,isAxisOn,nbins,plotGroup,~,extra] = internal.stats.parseArgs(paramNames,defaults,varargin{:});

if iscell(nbins)
    nbins = nbins{:};
end
if numel(nbins) == 1
   nbins = repmat(nbins,nGrp,1); 
end

if ~isempty(nbins) && numel(grp)~= numel(nbins)
    error(message('stats:internal:plotGroupedHist:BadBinNumber'));
end
if any(nbins<=0)
    error(message('stats:internal:plotGroupedHist:NegativeBinNumber'));
end

% calculate bin width
if isempty(nbins) && nGrp ~= 1
    for i = 1:nGrp
        xg = x(grpID == grp(i));
        [~,edges] = histcounts(xg,'BinMethod','scott','Norm','pdf');
        bw(i) = edges(2) - edges(1);
    end
    bw = min(bw);   % common bin width based on the group with the narrowest bins
    if range(x)/bw>20
        bw = range(x)/20; % limit the max number of bins(20)
    end
end

dataCursorBehaviorObj = hgbehaviorfactory('Datacursor');
set(dataCursorBehaviorObj,'UpdateFcn',@groupedhistDatatipCallback);

if nGrp == 1 
    if isempty(nbins)
        hXLines(1) = histogram(h,x,extra{:});
    elseif nGrp == 1 && ~isempty(nbins)
        hXLines(1) = histogram(h,x,nbins,extra{:});
    end
else
    for i = 1:nGrp
        xg = x(grpID == grp(i));
        if isempty(nbins)
           hXLines(i) = histogram(h,xg,'BinWidth',bw,extra{:});          
        else
           hXLines(i) = histogram(h,xg,nbins(i),extra{:});
        end
        hgaddbehavior(hXLines(i),dataCursorBehaviorObj);
        if ~isempty(gname)
            setappdata(hXLines(i),'groupname',gname{i});
        end
        hold on
    end
end
hold off

if isempty(clr)
    clr = lines(nGrp);
elseif ischar(clr) && isvector(clr)
    clr = clr(:);
end
if length(clr) == 3
    clr = [clr;clr];
end
clr = internal.stats.cycleLineProperties(nGrp,clr);

if plotGroup
    if  strcmpi(hXLines(1).DisplayStyle,'stairs')
        for i = 1:nGrp
            set(hXLines(i),'EdgeColor',clr(i,:));
        end
    else
        for i = 1:nGrp
            set(hXLines(i),'FaceColor',clr(i,:));
        end
    end
end

% Draw a black horizontal line as the X axis if axis is 'off'.
if ~isAxisOn
    line(h.XLim,[0 0],'Color','k','HandleVisibility','off');    
end

if nargout > 0
    hg = hXLines;
end
end

function datatipTxt = groupedhistDatatipCallback(~,evt)
target = get(evt,'Target');
pos = get(evt,'Position');
x = pos(1); 
y = pos(2);
bw = get(target,'Binwidth');
xL = x - bw/2;
xR = x + bw/2;
target = get(evt,'Target');
groupname = getappdata(target,'groupname');
datatipTxt = {['Value: ',num2str(y,4)]...
              ['BinEdges: [',num2str(xL),' ', num2str(xR),']']...
              ['Group: ',groupname]};
end




