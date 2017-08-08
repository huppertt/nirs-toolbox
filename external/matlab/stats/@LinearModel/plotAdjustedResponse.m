function hout=plotAdjustedResponse(model,var,varargin)
%plotAdjustedResponse Adjusted response plot for regression model
%   plotAdjustedResponse(LM,VAR) creates an adjusted response plot for the 
%   variable VAR in the LinearModel LM. The adjusted response plot shows
%   the fitted response as a function of VAR, with the other predictors
%   averaged out by averaging the fitted values over the data used in the
%   fit. Adjusted data points are computed by adding the residual to the
%   adjusted fitted value for each observation.
%
%   H = plotAdjustedResponse(...) returns a vector H of handles to the
%   lines in the plot.
%
%   The VAR argument can be followed by parameter/value pairs to specify
%   properties of the lines in the plot, such as 'LineWidth' and 'Color'.
%
%   The data cursor tool in the figure window will display the X and Y
%   values for any data point, along with the observation name or number.
%
%    Example:
%      % Model MPG as a function of Weight, and create a plot showing
%      % the effect of Weight averaged over the Year values
%      load carsmall
%      d = dataset(MPG,Weight);
%      d.Year = ordinal(Model_Year);
%      lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%      plotAdjustedResponse(lm,'Weight')
%
%   See also LinearModel, plotEffects, plotInteraction, plotAdded.

%   Copyright 2011-2013 The MathWorks, Inc.

% Plot adjusted response as function of predictor vnum
narginchk(2,Inf);
internal.stats.plotargchk(varargin{:});

% Get information about terms and predictors
terminfo = getTermInfo(model);
[xdata,vname,vnum] = getVar(model,var);

if ~model.Formula.InModel(vnum)
    if strcmp(vname,model.Formula.ResponseName)
        error(message('stats:LinearModel:ResponseNotAllowed', vname));
    else
        error(message('stats:LinearModel:NotPredictor', vname));
    end

elseif terminfo.isCatVar(vnum)
    % Compute fitted values for each level of this predictor
    [xdata,xlabels] = grp2idx(xdata);
    nlevels = length(xlabels);
    xi = (1:nlevels)';
else
    % Define a grid of values; combine with data
    xi = linspace(min(xdata),max(xdata))';
    nlevels = length(xi);
    xi = [xi; xdata];
end

% Compute adjusted fitted values as a function of this
% predictor
fxi = getAdjustedResponse(model,vnum,xi,terminfo);

% Separate values on grid from values at data points
if terminfo.isCatVar(vnum)
    d = double(xdata);
    fx(~isnan(d)) = fxi(d(~isnan(d)));
    fx(isnan(d)) = NaN;
    fx = fx(:);
else
    fx = fxi(nlevels+1:end);
    xi = xi(1:nlevels);
    fxi = fxi(1:nlevels);
end

% Plot the results
resid = model.Residuals.Raw;
h = plot(xdata,fx+resid,'ro', varargin{:});
ax = ancestor(h,'axes');
washold = ishold(ax);
if ~washold
    hold(ax,'on');
end
h = [h; plot(ax,xi,fxi,'b-')];
if ~washold
    hold(ax,'off');
end
set(h(1),'Tag','data');
set(h(2),'Tag','fit');
legend(ax,'Adjusted data','Adjusted fit','location','best');
xlabel(ax,vname);
ylabel(ax,sprintf('%s',getString(message('stats:LinearModel:sprintf_Adjusted',model.ResponseName))));
title(ax,getString(message('stats:LinearModel:title_AdjustedResponsePlot')))

if terminfo.isCatVar(vnum)
    set(ax,'XTick',xi,'XTickLabel',xlabels,'XLim',[.5,max(xi)+.5]);
end

% Define data tips
ObsNames = model.ObservationNames;
internal.stats.addLabeledDataTip(ObsNames,h(1),[]);

if nargout>0
    hout = h;
end
