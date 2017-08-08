function hout = plotInteraction(model,var1,var2,ptype)
%plotInteraction Plot interaction effects of two predictors.
%   plotInteraction(LM,VAR1,VAR2) creates a plot of the interaction effects
%   of the predictors VAR1 and VAR2 in the LinearModel LM. This plot shows
%   the estimated effect on the response from changing each predictor value
%   from one value to another with the effects of other predictors averaged
%   out. It also shows the estimated effect with the other predictor fixed
%   at certain values. The values are chosen to produce a relatively large
%   effect on the response. This plot enables you to examine whether the
%   effect of one predictor depends on the value of the other predictor.
%
%   Each effect is shown as a circle, with a horizontal bar showing the
%   confidence interval for the estimated effect. The effect values are
%   computed from the adjusted response curve, as shown by the
%   plotAdjustedResponse function.
%
%   plotInteraction(LM,VAR1,VAR2,PTYPE) with PTYPE='predictions' shows the
%   adjusted response curve as a function of VAR2, with VAR1 fixed at
%   certain values. The default form of the plot is produced with
%   PTYPE='effects'.
%
%   H = plotInteraction(...) returns a vector H of handles to the lines in
%   the plot.
%
%    Example:
%      % Plot the effect of one predictor in a regression model at
%      % fixed values of another predictor
%      load carsmall
%      d = dataset(MPG,Weight);
%      d.Year = ordinal(Model_Year);
%      lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%      subplot(1,2,1)
%      plotInteraction(lm,'Year','Weight','predictions')
%
%      % See how this changes if we add an interaction term
%      lm = addTerms(lm,'Weight*Year')
%      subplot(1,2,2)
%      plotInteraction(lm,'Year','Weight','predictions')
%
%   See also LinearModel, plotAdjustedResponse, plotEffects.

%   Copyright 2011-2013 The MathWorks, Inc.

%   For a 'predictions' plot, there is one value in H for each curve shown
%   on the plot. For an 'effects' plot, H(1) is the handle for the main
%   effects. H(2) is the handle for the confidence interval for the first
%   main effect, and H(3) is the confidence interval for the second main
%   effect. The remaining entries in H are handles for the conditional
%   effects and their confidence intervals. Handles associated with the
%   main effects have the tag 'main'. Handles associated with conditional
%   effects at the top and bottom have the tag 'conditional1' and
%   'conditional2', respectively.

if nargin<4
    ptype = 'effects';
end

% Plot interaction effects between two predictors

% Get information about terms and predictors
terminfo = getTermInfo(model);
[vname1,vnum1] = identifyVar(model,var1);
[vname2,vnum2] = identifyVar(model,var2);

if isequal(vname1,model.ResponseName) || isequal(vname2,model.ResponseName)
    error(message('stats:LinearModel:ResponseNotAllowed',model.ResponseName))
elseif isequal(vnum1,vnum2)
    error(message('stats:LinearModel:DifferentPredictors'))
end

switch(ptype)
    case 'effects'
        h = plotInteractionEffects(model,vnum1,vnum2,vname1,vname2,terminfo);
    case 'predictions'
        subset = model.ObservationInfo.Subset;
        xdata1 = model.Variables.(vname1);
        xdata1 = xdata1(subset,:);
        xdata2 = model.Variables.(vname2);
        xdata2 = xdata2(subset,:);
        h = plotInteractionPredictions(model,vnum1,vnum2,vname1,vname2,xdata1,xdata2,terminfo);
    otherwise
        error(message('stats:LinearModel:BadEffectsType'));
end

if nargout>0
    hout = h;
end

% ------------------
function h = plotInteractionPredictions(model,vnum1,vnum2,vname1,vname2,xdata1,xdata2,terminfo)

if terminfo.isCatVar(vnum1)
    % Compute fitted values for each level of this predictor
    [~,xlabels1] = grp2idx(xdata1);
    ngrid1 = length(xlabels1);
    xi1 = (1:ngrid1)';
else
    % Define a grid of values
    ngrid1 = 3;
    xi1 = linspace(min(xdata1),max(xdata1),ngrid1)';
    xlabels1 = cellstr(strjust(num2str(xi1),'left'));
end

if terminfo.isCatVar(vnum2)
    % Compute fitted values for each level of this predictor
    [~,xlabels2] = grp2idx(xdata2);
    ngrid2 = length(xlabels2);
    xi2 = (1:ngrid2)';
    plotspec = '-o';
else
    % Define a grid of values
    ngrid2 = 100;
    xi2 = linspace(min(xdata2),max(xdata2),ngrid2)';
    plotspec = '-';
end

y = zeros(ngrid2,length(xi1));
xi = [xi2 xi2];
for j = 1:size(y,2)
    xi(:,1) = xi1(j);
    y(:,j) = getAdjustedResponse(model,[vnum1 vnum2],xi,terminfo);
end

h = plot(1,1, xi2,y,plotspec, 'LineWidth',2);
set(h(1),'LineStyle','none','Marker','none','XData',[],'YData',[]);
title(sprintf('%s',getString(message('stats:LinearModel:sprintf_InteractionOfAnd',vname1,vname2))));
xlabel(vname2);
ylabel(sprintf('%s',getString(message('stats:LinearModel:sprintf_Adjusted',model.ResponseName))));
legend(h,[vname1; xlabels1(:)]);

if terminfo.isCatVar(vnum2)
    set(gca,'XTick',1:ngrid2,'XTickLabel',xlabels2,'XLim',[0.5,ngrid2+0.5]);
end
h(1) = [];  % special line used only to create legend label

% ------------------
function h = plotInteractionEffects(model,vnum1,vnum2,vname1,vname2,terminfo)
[effect,effectSE,effectName,x] = getEffects(model,[vnum1 vnum2],terminfo);
ci = [effect effect] + effectSE*tinv([.025 .975],model.DFE);

[ceffect1,ceffect1SE,ceffect1Name] = getConditionalEffect(model,vnum1,vnum2,x(1,:)',terminfo);
ci1 = [ceffect1 ceffect1] + ceffect1SE*tinv([.025 .975],model.DFE);

[ceffect2,ceffect2SE,ceffect2Name] = getConditionalEffect(model,vnum2,vnum1,x(2,:)',terminfo);
ci2 = [ceffect2 ceffect2] + ceffect2SE*tinv([.025 .975],model.DFE);

% Plot the results
gap = 2;
y0 = [1; 2+gap+length(ceffect1)];
y1 = 1 + (1:length(ceffect1))';
y2 = 2 + gap + length(ceffect1) + (1:length(ceffect2))';
y = [y0(1); y1; y0(2); y2];
allnames = [effectName(1); ceffect1Name; effectName(2);ceffect2Name];

h = plot(effect,y0,'bo', ci',[y0 y0]','b-', 'LineWidth',2,'Tag','main');
washold = ishold;
hold on
h = [h;
     plot(ceffect1,y1,'Color','r','LineStyle','none','Marker','o','Tag','conditional1'); ...
     plot(ci1',[y1 y1]','Color','r','Tag','conditional1'); ...
     plot(ceffect2,y2,'Color','r','LineStyle','none','Marker','o','Tag','conditional2'); ...
     plot(ci2',[y2 y2]','Color','r','Tag','conditional2')];
if ~washold
    hold off
end
title(sprintf('%s',getString(message('stats:LinearModel:sprintf_InteractionOfAnd',vname1,vname2))));
xlabel(getString(message('stats:LinearModel:xylabel_Effect')));
set(gca,'YTick',y,'YTickLabel',allnames,'YLim',[.5,max(y)+.5],'YDir','reverse');
dfswitchyard('vline',gca,0,'LineStyle',':','Color','k');

