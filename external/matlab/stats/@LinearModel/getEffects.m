function [effects,effectSEs,effectnames,effectXs] = getEffects(model,vars,terminfo)
% Get the main effect of each specified predictor, computed as the maximum
% change in adjusted response between two different predictor values.

%   Copyright 2011 The MathWorks, Inc.

if nargin<3
    % Get information about terms and predictors
    terminfo = getTermInfo(model);
end
if nargin<2
    vars = model.PredictorNames;
end

npred = length(vars);
effectnames = cell(npred,1);
effects = zeros(npred,1);
effectSEs = zeros(npred,1);
effectXs = zeros(npred,2);

for j=1:length(vars)
    [xdata,vname,vnum] = getVar(model,vars(j));
    
    if terminfo.isCatVar(vnum)
        % Compute fitted values for each level of this predictor
        [~,xlabels] = grp2idx(xdata);
        xi = (1:length(xlabels))';
    else
        % Define a grid of values
        xi = linspace(min(xdata),max(xdata),101)';
    end
    
    % Compute adjusted fitted values as a function of this predictor
    [fxi,fxiVar] = getAdjustedResponse(model,vnum,xi,terminfo);
    
    % Compute main effect, its standard error, and its label
    [maxf,maxloc] = max(fxi);
    [minf,minloc] = min(fxi);
    effect = maxf - minf;
    effectSE = sqrt(max(0,fxiVar(minloc,minloc) + fxiVar(maxloc,maxloc) - 2*fxiVar(minloc,maxloc)));
    if terminfo.isCatVar(vnum)
        effectname = sprintf('%s',getString(message('stats:LinearModel:sprintf_EffectAtoB',vname,xlabels{minloc},xlabels{maxloc})));
    else
        if minloc>maxloc
            effect = -effect;
            temp = minloc;
            minloc = maxloc;
            maxloc = temp;
        end
        effectname = sprintf('%s',getString(message('stats:LinearModel:sprintf_EffectAtoB',vname,num2str(xi(minloc)),num2str(xi(maxloc)))));
    end
    
    effectX = [xi(minloc), xi(maxloc)];
    
    effects(j) = effect;
    effectnames{j} = effectname;
    effectSEs(j) = effectSE;
    effectXs(j,:) = effectX;
end
