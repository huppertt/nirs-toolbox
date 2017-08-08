function hout = plot(lm,varargin)
%PLOT Summary plot of regression model
%   PLOT(LM) produces a summary plot of the LinearModel LM.
%
%   If LM has exactly one predictor, the plot is a plot of the response as
%   a function of the predictor, with the fit and confidence bounds
%   superimposed.
%
%   If LM has more than one predictor, the plot is an added variable plot
%   for the model as a whole. If LM has no predictors, the plot is a
%   histogram of residuals.
%
%   H = PLOT(LM) returns a vector H of handles to the lines in the plot.
%
%   The LM argument can be followed by parameter/value pairs to specify
%   properties of the lines in the plot, such as 'LineWidth' and 'Color'.
%
%    Example:
%       % Plot of census data with fitted quadratic curve
%       load census
%       lm = fitlm(cdate,pop,'quadratic');
%       plot(lm)
%
%   See also LinearModel, plotAdded.

%   Copyright 2011-2013 The MathWorks, Inc.

p = length(lm.PredictorNames);
internal.stats.plotargchk(varargin{:});

if p==0
    % No predictors, plot residuals
    h = plotResiduals(lm,'histogram');
elseif p==1
    % One predictor, plot x vs y with fit and confidence bounds
    h = plotxy(lm,varargin{:});
else
    h = plotAdded(lm,[],varargin{:});
end

if nargout>0
    hout = h;
end

function h=plotxy(lm,varargin)

% Only one predictor, where is it?
col = lm.PredLocs;
xname = lm.PredictorNames{1};

% Get its values
xdata = getVar(lm,col);
y = getResponse(lm);
ObsNames = lm.ObservationNames;

iscat = lm.VariableInfo.IsCategorical(col);

if iscat
    % Compute fitted values for each level of this predictor
    [x,xlabels,levels] = grp2idx(xdata);
    tickloc = (1:length(xlabels))';
    ticklab = xlabels;
    xx = tickloc;
else
    x = xdata;
    xx = linspace(min(x), max(x))';
    levels = xx;
end
nlevels = size(levels,1);

% Make sure NaNs match up to avoid having unused values (those paired with
% NaN in the other variable) affect the plot.
t = isnan(x) | isnan(y);
if any(t)
    x(t) = NaN;
    y(t) = NaN;
end

% Predict and plot
if isa(lm.Variables,'dataset') || isa(lm.Variables,'table')
    % Create table to hold this variable
    X = lm.Variables(ones(nlevels,1),:);
    X.(xname) = levels(:);       % for prediction
else
    % Create matrix to hold this variable
    npreds = lm.NumVariables-1;
    X = zeros(length(xx),npreds);
    X(:,col) = xx;
end 
[yfit,yci] = lm.predict(X);
h = plot(x,y,'bx', varargin{:});
ax = ancestor(h,'axes');
washold = ishold(ax);
hold(ax,'on')
h = [h; plot(ax,xx,yfit,'r-' ,xx,yci,'r:')];
if ~washold
    hold(ax,'off')
end

if iscat
    set(ax,'XTick',tickloc','XTickLabel',ticklab);
    set(ax,'XLim',[tickloc(1)-0.5, tickloc(end)+0.5]);
end

yname = lm.ResponseName;
title(ax,sprintf('%s',getString(message('stats:LinearModel:sprintf_AvsB',yname,xname))),'Interpreter','none');
set(xlabel(ax,xname),'Interpreter','none');
set(ylabel(ax,yname),'Interpreter','none');
legend(ax,h(1:3),getString(message('stats:LinearModel:legend_Data')), ...
                 getString(message('stats:LinearModel:legend_Fit')), ...
                 getString(message('stats:LinearModel:legend_ConfidenceBounds')), ...
                 'location','best')

% Define data tips
internal.stats.addLabeledDataTip(ObsNames,h(1),h(2:end));
