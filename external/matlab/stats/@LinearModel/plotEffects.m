function hout = plotEffects(model)
%plotEffects Plot main effects of each predictor.
%   plotEffects(LM) produces an effects plot for the predictors in the
%   LinearModel LM. This plot shows the estimated effect on the response
%   from changing each predictor value from one value to another. The two
%   values are chosen to produce a relatively large effect on the response.
%   This plot enables you to compare the effects of the predictors,
%   regardless of the scale of the predictor measurements.
%
%   Each effect is shown as a circle, with a horizontal bar showing the
%   confidence interval for the estimated effect. The effect values are
%   computed from the adjusted response curve, as shown by the
%   plotAdjustedResponse function.
%
%   H = plotEffects(LM) returns a vector H of handles to the lines in the
%   plot. H(1) is a handle to the line that defines the effect estimates,
%   and H(J+1) is a handle to the line that defines the confidence interval
%   for the effect of predictor J.
%
%    Example:
%      % Plot the effects of two predictors in a regression model
%      load carsmall
%      d = dataset(MPG,Weight);
%      d.Year = ordinal(Model_Year);
%      lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%      plotEffects(lm)
%
%      % Verify that plotted effect of Weight matches what we would
%      % calculate by evaluating the fitted model
%      feval(lm,4732,'70') - feval(lm,1795,'70')
%
%   See also LinearModel, plotAdjustedResponse, plotInteraction.

%   Copyright 2011-2013 The MathWorks, Inc.

% Plot main effects of each predictor
[effect,effectSE,effectname] = getEffects(model);

% Plot the results
y = (1:length(effect))';
ci = [effect effect] + effectSE*tinv([.025 .975],model.DFE);
h = plot(effect,y,'bo', ci',[y y]','b-');
set(h(1),'Tag','estimate');
set(h(2:end),'Tag','ci');
xlabel(getString(message('stats:LinearModel:xylabel_MainEffect')));
set(gca,'YTick',y,'YTickLabel',effectname,'YLim',[.5,max(y)+.5],'YDir','reverse');
dfswitchyard('vline',gca,0,'LineStyle',':','Color','k');

if nargout>0
    hout = h;
end
