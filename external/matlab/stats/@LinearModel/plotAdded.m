function hout = plotAdded(model,cnum,varargin)
%plotAdded Added variable plot
%    plotAdded(LM,CNUM) produces an added variable plot for the term that
%    multiples coefficient number CNUM in the LinearModel LM. An added
%    variable plot illustrates the incremental effect on the response of
%    this term by removing the effects of all other terms. The slope of the
%    fitted line is coefficient number CNUM.
%
%    plotAdded(LM,V) where V is a vector of coefficient numbers generalizes
%    the added variable plot to multiple terms. It plots the response vs.
%    the best linear combination of the specified terms, after adjusting
%    for all other terms. The slope of the fitted line is the coefficient
%    of the linear combination of the specified terms projected onto the
%    best-fitting direction.
%
%    plotAdded(LM) produces a generalized added variable plot for all terms
%    except the constant term.
%
%    The data cursor tool in the figure window will display the X and Y
%    values for any data point, along with the observation name or number.
%
%    Example:
%      % Model MPG as a function of Weight, and create an added variable
%      % plot to view the marginal effect of the squared term
%      load carsmall
%      d = dataset(MPG,Weight);
%      d.Year = ordinal(Model_Year);
%      lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%      plotAdded(lm,'Weight^2')
%
%    See also LinearModel, addedvarplot.

%   Copyright 2011-2013 The MathWorks, Inc.

%    Reference:
%    Sall, J. (1990), "Leverage plots for general linear hypotheses,"
%    American Statistician, v. 44, pp. 308-315. The plot produced by this
%    function is a modified version of Sall's plot. Here the x-axis
%    variable is scaled so that the fitted line has slope given by the
%    coefficients rather than slope 1. This variant is mentioned in
%    section 2.2 of the paper for the case of a single predictor.

internal.stats.plotargchk(varargin{:});

% Get model info in a form expected by addedvarplot
sub = model.ObservationInfo.Subset;
stats.source = 'stepwisefit';
stats.B = model.Coefs;
stats.SE = model.CoefSE;
stats.dfe = model.DFE;
stats.covb = model.CoefficientCovariance;
% need this for "out" coefficients xr = stats.xr(:,outnum);
stats.yr = model.Residuals.Raw(sub);
stats.wasnan = ~sub;
stats.wts = get_CombinedWeights_r(model);
stats.mse = model.MSE;
[~,p] = size(model.Design);

% Check or create CNUM
terminfo = getTermInfo(model);
constrow = find(all(terminfo.terms==0,2),1);
if isempty(constrow)
    constrow = NaN;
end
ncoefs = length(model.Coefs);
if nargin<2 || isempty(cnum)
    % Find the number for all coefficients except the intercept
    cnum = find(terminfo.designTerms~=constrow);
end
if isrow(cnum) && ischar(cnum)
    termnum =  find(strcmp(model.Formula.TermNames,cnum));
    if isscalar(termnum)
        cnum = find(terminfo.designTerms==termnum);
    else
        cnum = find(strcmp(model.CoefficientNames,cnum));
        if ~isscalar(cnum)
            error(message('stats:LinearModel:BadCoefficientName'));
        end
    end
elseif isempty(cnum) || ~isvector(cnum) || ~all(ismember(cnum,1:ncoefs))
    error(message('stats:LinearModel:BadCoefficientNumber'));
end
cnum = sort(cnum);
if ~isscalar(cnum) && any(diff(cnum)==0)
    error(message('stats:LinearModel:RepeatedCoeffients'));
end
   
% Create added variable plot
y = getResponse(model);
h = addedvarplot(model.Design(sub,:),y(sub),cnum,true(1,p),stats,[],false,varargin{:});

% Customize axis properties
ax = ancestor(h(1),'axes');
ylabel(ax,sprintf('%s',getString(message('stats:LinearModel:sprintf_Adjusted',model.ResponseName))),'Interpreter','none');
tcols = terminfo.designTerms(cnum);
if isscalar(cnum) % single coefficient
    thetitle = sprintf('%s',getString(message('stats:LinearModel:title_AddedVariablePlotFor',model.CoefficientNames{cnum})));
    thexlabel = sprintf('%s',getString(message('stats:LinearModel:sprintf_Adjusted',model.CoefficientNames{cnum})));
elseif ~any(tcols==constrow) && length(cnum)==ncoefs-1 % all except const
    thetitle = getString(message('stats:LinearModel:title_AddedVariablePlotModel'));
    thexlabel = getString(message('stats:LinearModel:xylabel_AdjustedWholeModel'));
elseif all(tcols==tcols(1)) && length(tcols)==sum(terminfo.designTerms==tcols(1))
    % all coefficients for one term
    thetitle = sprintf('%s',getString(message('stats:LinearModel:title_AddedVariablePlotFor',model.Formula.TermNames{tcols(1)})));
    thexlabel = sprintf('%s',getString(message('stats:LinearModel:sprintf_Adjusted',model.Formula.TermNames{tcols(1)})));
else
    thetitle = getString(message('stats:LinearModel:title_AddedVariablePlotTerms'));
    thexlabel = getString(message('stats:LinearModel:xylabel_AdjustedSpecifiedTerms'));
end
title(ax,thetitle,'Interpreter','none');
xlabel(ax,thexlabel,'Interpreter','none');

% Define data tips
ObsNames = model.ObservationNames;
internal.stats.addLabeledDataTip(ObsNames,h(1),h(2:end));

if nargout>0
    hout = h;
end
