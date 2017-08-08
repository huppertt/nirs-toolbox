function [effect,effectSE,effectName] = getConditionalEffect(model,var1,var2,xi1,terminfo)
% Get effect of v1 conditional on various values of v2. xi1 specifies two
% values of v1, with the effect defined as the difference in adjusted
% response between those values. Typically xi1 comes from getEffects.

%   Copyright 2011 The MathWorks, Inc.

if nargin<5
    % Get information about terms and predictors
    terminfo = getTermInfo(model);
end
[~,vnum1] = identifyVar(model,var1);
[xdata2,vname2,vnum2] = getVar(model,var2);

if terminfo.isCatVar(vnum2)
    % Compute fitted values for each level of this predictor
    [~,xlabels2] = grp2idx(xdata2);
    ngrid = length(xlabels2);
    xi2 = (1:ngrid)';
else
    % Define a grid of values
    ngrid = 3;
    xi2 = linspace(min(xdata2),max(xdata2),ngrid)';
end

xi = [repmat(xi1(1),ngrid,1) xi2; ...
      repmat(xi1(2),ngrid,1) xi2];

% Compute adjusted fitted values as a function of this
% predictor
[fxi,fxiVar] = getAdjustedResponse(model,[vnum1 vnum2],xi,terminfo);

% Compute conditional effect for predictor 1, its standard error, and its label
effect = fxi(ngrid+1:2*ngrid) - fxi(1:ngrid);
fxiVarDiag = diag(fxiVar);
fxiCov = diag(fxiVar(1:ngrid,ngrid+1:2*ngrid));
effectSE = sqrt(max(fxiVarDiag(1:ngrid) + fxiVarDiag(ngrid+1:2*ngrid) - 2*fxiCov,0));
if terminfo.isCatVar(vnum2)
    effectName = strcat(sprintf('%s=',vname2),xlabels2(:));
else
    effectName = {sprintf('%s=%g',vname2,xi2(1)); ...
                  sprintf('%s=%g',vname2,xi2(2)); ...
                  sprintf('%s=%g',vname2,xi2(3))};
end
