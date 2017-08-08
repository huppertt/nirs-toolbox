function hout = plotResiduals(model,plottype,varargin)
%plotResiduals Plot residuals of fitted model
%    plotResiduals(M,PLOTTYPE) plots the residuals from model M in
%    a plot of type PLOTTYPE. Valid values for PLOTTYPE are:
%
%       'caseorder'     residuals vs. case (row) order
%       'fitted'        residuals vs. fitted values
%       'histogram'     histogram (default)
%       'lagged'        residual vs. lagged residual (r(t) vs. r(t-1))
%       'probability'   normal probability plot
%       'symmetry'      symmetry plot
%
%    H = plotResiduals(...) returns a handle to the lines or patches in the
%    plot.
%
%    The PLOTTYPE argument can be followed by parameter/value pairs to
%    specify additional properties of the primary line in the plot. For
%    example, plotResiduals(M,'fitted','Marker','s') uses a square marker.
%
%    For many of these plots, the data cursor tool in the figure window
%    will display the X and Y values for any data point, along with the
%    observation name or number.
%
%    See also classreg.regr.modelutils.plotDiagnostics.

%   Copyright 2011-2014 The MathWorks, Inc.

% Plot residuals in any number of ways
narginchk(1,Inf);
internal.stats.plotargchk(varargin{:});

% Get plot type
if nargin<2
    plottype = 'histogram';
end

alltypes = {'caseorder' 'fitted' 'histogram' 'lagged' 'probability' 'symmetry'};
tf = strncmpi(plottype,alltypes,length(plottype));
if sum(tf)~=1
    error(message('stats:LinearModel:BadResidualPlotType'));
end
plottype = alltypes{tf};

if isempty(varargin)
    residtype = 'Raw';
    % wantWeighted is true if you want weighted residuals.
    wantWeighted = false;
else
    [residtype,wantWeighted,~,args] = internal.stats.parseArgs({'residualtype','weighted'},{'Raw',false},varargin{:});
    varargin = args;
end

% Get residuals and observation names from the fit
if ( wantWeighted )
    r = model.WeightedResiduals.(residtype);
else
    r = model.Residuals.(residtype);
end
ObsNames = model.ObservationNames;

h = [];
h0 = [];
labels = [];

switch(plottype)
    
    case 'caseorder'
        h = plot(r,'bx',varargin{:});
        ax = ancestor(h,'axes');
        xlim = get(ax,'XLim');
        h0 = line(xlim,[0 0],'Color','k','LineStyle',':','XLimInclude','off','Parent',ax);
        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_RowNumber')))
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_Residuals')))
        title(ax,getString(message('stats:classreg:regr:modelutils:title_CaseOrderPlotOfResiduals')))

    case 'histogram'
        h = internal.stats.histogram(r,'rule','scott','type','Probability',varargin{:});
        ax = ancestor(h,'axes');
        title(ax,getString(message('stats:classreg:regr:modelutils:title_HistogramOfResiduals')))
        
    case 'fitted'
        yhat = predict(model);
        yhat(isnan(r)) = NaN;
        h = plot(yhat,r,'bx',varargin{:});
        ax = ancestor(h,'axes');
        xlim = get(ax,'XLim');
        h0 = line(xlim,[0 0],'Color','k','LineStyle',':','XLimInclude','off','Parent',ax);
        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_FittedValues')))
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_Residuals')))
        title(ax,getString(message('stats:classreg:regr:modelutils:title_PlotOfResidualsVsFittedValues')))
        
    case 'lagged'
        h = plot(r(1:end-1),r(2:end),'bx',varargin{:});
        ax = ancestor(h,'axes');
        
        xlim = get(ax,'XLim');
        h0 = line(xlim,[0 0],'Color','k','LineStyle',':','XLimInclude','off','Parent',ax);
        
        ylim = get(ax,'YLim');
        h0(2) = line([0 0],ylim,'Color','k','LineStyle',':','YLimInclude','off','Parent',ax);
        
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_Residualt')))
        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_Residualt1')))
        title(ax,getString(message('stats:classreg:regr:modelutils:title_PlotOfResidualsVsLaggedResiduals')))
        labels = @(j) laggedDataTip(j);

    case 'probability'
        [ax,args] = axescheck(varargin{:});
        if isempty(ax)
            h = probplot(r);
        else
            h = probplot(ax,r);
        end
        ax = ancestor(h(1),'axes');
        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_Residuals')))
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_Probability')))
        title(ax,getString(message('stats:classreg:regr:modelutils:title_NormalProbabilityPlotOfResiduals')))
        inds = getappdata(h(1),'originds');
        setappdata(h(1),'originds',[]); % mark as not from a dataset
        if isempty(ObsNames)
            labels = strcat({'Observation '},strjust(num2str(inds(:)),'left'));
        else
            labels = ObsNames(inds);
        end
        if ~isempty(args)
            set(h(1),args{:});
        end
        h0 = h(2:end);
        h = h(1);
        
    case 'symmetry'
        [sr,idx] = sort(r);
        n = sum(~isnan(sr));
        med = median(sr(1:n));
        nlo = floor(n/2);
        nhi = n+1-nlo;
        h = plot(med-sr(nlo:-1:1), sr(nhi:n)-med, 'bx',varargin{:});
        ax = ancestor(h,'axes');
        
        xlim = get(ax,'XLim');
        ylim = get(ax,'YLim');
        xymax = max(xlim(2),ylim(2));
        h0 = line([0 xymax],[0 xymax],'Color','k','LineStyle',':',...
                                      'YLimInclude','off','XLimInclude','off','Parent',ax);
        
        xlabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_LowerTail')))
        ylabel(ax,getString(message('stats:classreg:regr:modelutils:xylabel_UpperTail')))
        title(ax,getString(message('stats:classreg:regr:modelutils:title_SymmetryPlotOfResidualsAroundTheirMedian')))
        labels = @(j) symmetryDataTip(j,idx,nlo,nhi);
end

if ~isempty(h)
    set(h,'Tag','main');
end
if ~isempty(h0)
    set(h0,'Tag','reference');
end
% Define data tips
if ~isempty(h) && ~strcmp(plottype,'histogram')
    if isempty(labels)
        labels = ObsNames;
    end
    internal.stats.addLabeledDataTip(labels,h,h0);
end

if nargout>0
    hout = [h(:); h0(:)];
end

% -----------
function txt = laggedDataTip(j)
txt = getString(message('stats:classreg:regr:modelutils:sprintf_XvsY',getObsName(j),getObsName(j+1)));
end

function txt = symmetryDataTip(j,idx,nlo,nhi)
txt = getString(message('stats:classreg:regr:modelutils:sprintf_XvsY',getObsName(idx(nlo+1-j)),getObsName(idx(nhi+j-1))));
end

function obsname = getObsName(j)
if iscell(ObsNames) && ~isempty(ObsNames)
    obsname = ObsNames{j};
else
    obsname = sprintf('%s',getString(message('stats:classreg:regr:modelutils:sprintf_ObservationNumber',j)));
end
end

end
