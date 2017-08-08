function hout = plotDiagnostics(model,plottype,varargin)
%plotDiagnostics Plot diagnostics of fitted model
%    plotDiagnostics(M,PLOTTYPE) plots diagnostics from model M in
%    a plot of type PLOTTYPE. The default value of PLOTTYPE is 'leverage'.
%    Valid values for PLOTTYPE are:
%
%       'contour'      residual vs. leverage with overlaid Cook's contours
%       'cookd'        Cook's distance
%       'covratio'     delete-1 ratio of determinant of covariance
%       'dfbetas'      scaled delete-1 coefficient estimates
%       'dffits'       scaled delete-1 fitted values
%       'leverage'     leverage (diagonal of Hat matrix)
%       's2_i'         delete-1 variance estimate
%
%    H = plotDiagnostics(...) returns handles to the lines in the plot.
%
%    The PLOTTYPE argument can be followed by parameter/value pairs to
%    specify additional properties of the primary line in the plot. For
%    example, plotDiagnostics(M,'leverage','Marker','s') uses a square
%    marker.
%
%    The data cursor tool in the figure window will display the X and Y
%    values for any data point, along with the observation name or number.
%    It also displays the coefficient name for 'dfbetas'.
%
%    See also classreg.regr.modelutils.plotResiduals.

%   Copyright 2011-2015 The MathWorks, Inc.

% Plot diagnostics in any number of ways
narginchk(1,Inf);

% Get plot type
if nargin<2
    plottype = 'leverage';
end
internal.stats.plotargchk(varargin{:});

alltypes = {'contour' 'cookd' 'covratio' 'dfbetas' 'dffits' 'leverage' 's2_i'};
tf = strncmpi(plottype,alltypes,length(plottype));
if sum(tf)~=1
    error(message('stats:LinearModel:BadDiagnosticsPlotType'));
end
plottype = alltypes{tf};

% Prepare default values
h1 = [];           % handles to receive data tips
h0 = [];           % handles not to receive data tips
linelabels = [];   % labels for data tip of each plotted line in H1

n = model.NumObservations;
p = model.NumEstimatedCoefficients;

switch(plottype)
    
    case 'contour'  % Contour plot of components of Cook's distance
        % Scatter plot of residuals vs. leverage
        x = model.Diagnostics.Leverage;
        y = model.Residuals.Raw;
        h1 = plot(x,y,'rx','LineWidth',2,varargin{:},'DisplayName','Observations');
        ax = ancestor(h1,'axes');
        
        % Cook's distance is a function of (x,y), so add contours
        xlim = get(ax,'XLim');
        ylim = get(ax,'YLim');
        sigma = varianceParam(model);
        X = linspace(max(.01,xlim(1)),xlim(2),31);
        Y = linspace(ylim(1),ylim(2),30);
        Z = bsxfun(@times,abs(Y)'.^2, (X./(1-X).^2)) / (p*sigma^2);
        washold = ishold(ax);
        if ~washold
            hold(ax,'on');
        end
        grey = .7 * [1 1 1];
        v = getNiceContours(model.Diagnostics.CooksDistance);
        [C,h0] = contour(ax,X,Y,Z,v,'LineStyle',':','Color',grey);
        clabel(C,h0,'Color',grey);
        set(h0,'DisplayName',getString(message('stats:classreg:regr:modelutils:CooksDistanceContours')));
        if ~washold
            hold(ax,'off');
        end
        
        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_Leverage')))
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_Residual')))
        title(ax,getString(message('stats:classreg:regr:modelutils:title_CooksDistanceFactorization')))

    case 'cookd'  % Plot of Cook's distance vs row number
        y = model.Diagnostics.CooksDistance;
        subset = model.ObservationInfo.Subset;
        h1 = plot(y,'rx','DisplayName',getString(message('stats:classreg:regr:modelutils:xylabel_CooksDistance')),varargin{:});
        ax = ancestor(h1,'axes');
        xlim = get(ax,'XLim');
        
        % Reference line may help to judge influential points
        yref = 3*mean(y(isfinite(y) & subset));
        h0 = line(xlim,[yref yref],'Color','k','LineStyle',':','XLimInclude','off','Parent',ax,'DisplayName',getString(message('stats:classreg:regr:modelutils:ReferenceLine')));
        
        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_RowNumber')))
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_CooksDistance')))
        title(getString(message('stats:classreg:regr:modelutils:title_CaseOrderPlotOfCooksDistance')))
        
    case 'covratio'  % Plot of covariance ratio vs row number
        y = model.Diagnostics.CovRatio;
        h1 = plot(y,'rx','DisplayName',getString(message('stats:classreg:regr:modelutils:xylabel_CovarianceRatio')),varargin{:});
        ax = ancestor(h1,'axes');
        xlim = get(ax,'XLim');

        % Reference lines may help to judge influential points
        yref = 1 + [-1 1]*3*p/n;
        h0 = line([xlim';NaN;xlim'],[yref([1 1])';NaN;yref([2 2])'],'Color','k','LineStyle',':','XLimInclude','off','Parent',ax,'DisplayName',getString(message('stats:classreg:regr:modelutils:ReferenceLine')));

        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_RowNumber')))
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_CovarianceRatio')))
        title(ax,getString(message('stats:classreg:regr:modelutils:title_CaseOrderPlotOfCovarianceRatio')))
        
    case 'dfbetas'  % Plot of scaled change in coefficients vs row number
        y = model.Diagnostics.Dfbetas;
        h1 = plot(y,'x',varargin{:});
        set(h1,{'DisplayName'},model.CoefficientNames');
        ax = ancestor(h1(1),'axes');
        xlim = get(ax,'XLim');

        % Reference lines may help to judge influential points
        yref = [-1 1]*3/sqrt(n);
        h0 = line([xlim';NaN;xlim'],[yref([1 1])';NaN;yref([2 2])'],'Color','k','LineStyle',':','XLimInclude','off','Parent',ax,'DisplayName',getString(message('stats:classreg:regr:modelutils:ReferenceLine')));
        
        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_RowNumber')))
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_ScaledChangeInCoefficients')))
        title(ax,getString(message('stats:classreg:regr:modelutils:title_CaseOrderPlotOfScaledChangeInCoefficients')))
        linelabels = model.CoefficientNames;

    case 'dffits'  % Plot of scaled change in fitted values vs row number
        y = model.Diagnostics.Dffits;
        h1 = plot(y,'rx','DisplayName',getString(message('stats:classreg:regr:modelutils:xylabel_ScaledChangeInFit')),varargin{:});
        ax = ancestor(h1,'axes');
        xlim = get(ax,'XLim');

        % Reference lines may help to judge influential points
        yref = [-1 1]*2*sqrt(p/n);
        h0 = line([xlim';NaN;xlim'],[yref([1 1])';NaN;yref([2 2])'],'Color','k','LineStyle',':','XLimInclude','off','Parent',ax,'DisplayName',getString(message('stats:classreg:regr:modelutils:ReferenceLine')));

        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_RowNumber')))
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_ScaledChangeInFit')))
        title(ax,getString(message('stats:classreg:regr:modelutils:title_CaseOrderPlotOfScaledChangeInFit')))
        
    case 'leverage'  % Plot of leverage (hat matrix diagonal) vs row number
        y = model.Diagnostics.Leverage;
        h1 = plot(y,'rx','DisplayName',getString(message('stats:classreg:regr:modelutils:xylabel_Leverage')),varargin{:});
        ax = ancestor(h1,'axes');
        xlim = get(ax,'XLim');
        
        % Reference line may help to judge influential points
        yref = 2*p/n;
        h0 = line(xlim,[yref yref],'Color','k','LineStyle',':','XLimInclude','off','Parent',ax,'DisplayName',getString(message('stats:classreg:regr:modelutils:ReferenceLine')));
        
        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_RowNumber')))
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_Leverage')))
        title(ax,getString(message('stats:classreg:regr:modelutils:title_CaseOrderPlotOfLeverage')))

    case 's2_i'  % Plot of leave-one-out variance vs row number
        y = model.Diagnostics.S2_i;
        h1 = plot(y,'rx','DisplayName',getString(message('stats:classreg:regr:modelutils:xylabel_LeaveoneoutVariance')),varargin{:});
        ax = ancestor(h1,'axes');
        xlim = get(ax,'XLim');

        % Reference line is the model MSE; no cutoff for extreme values
        yref = model.MSE;
        h0 = line(xlim,[yref yref],'Color','k','LineStyle',':','XLimInclude','off','Parent',ax,'DisplayName',getString(message('stats:classreg:regr:modelutils:MSE')));
        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_RowNumber')))
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_LeaveoneoutVariance')))
        title(ax,getString(message('stats:classreg:regr:modelutils:title_CaseOrderPlotOfLeaveoneoutVariance')))
end

% Define data tips
if ~isempty(h1)
    ObsNames = model.ObservationNames;
    internal.stats.addLabeledDataTip(ObsNames,h1,h0,linelabels);
end

if nargout>0
    hout = [h1(:); h0(:)];
end

% ------------------
function v = getNiceContours(vec)

% vec is a vector of non-negative values
vec = vec(isfinite(vec));
maxval = max(vec);

powOfTen = 10^floor(log10(maxval)); % next lower power of 10
relSize = maxval / powOfTen; % guaranteed in [1, 10)
if  relSize < 1.5
    v = powOfTen * (1:6)/4;
elseif relSize < 2.5
    v = powOfTen * (1:5)/2;
elseif relSize < 4
    v = powOfTen * (1:8)/2;
elseif relSize < 7.5
    v = powOfTen * (1:8);
else
    v = powOfTen * (1:5)*2;
end


