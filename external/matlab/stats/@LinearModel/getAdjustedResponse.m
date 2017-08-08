function [fxi,fxiVar] = getAdjustedResponse(model,var,xi,terminfo)
% Compute adjusted response as a function of a predictor

%   Copyright 2011 The MathWorks, Inc.

if nargin<4
    % Get information about terms and predictors
    terminfo = getTermInfo(model);
end
if isnumeric(var) % may be scalar or vector
    vnum = var;
else
    [~,vnum] = identifyVar(model,var);
end

% Remove this variable from all terms, and get the mean of the
% remaining term, plus a vector indicating this variable's power
% (continous) or setting (categorical)
[xrow,psmatrix,psflag] = reduceterm(model,vnum,terminfo);

% Create a matrix defining linear combinations of the coefficients
nrows = size(xi,1);
X = repmat(xrow,nrows,1);
for k=1:length(psflag)
    if psflag(k)
        % One row for each level of a categorical predictor; loop over levels
        for j=1:max(psmatrix(k,:))
            t = (psmatrix(k,:)==j);
            if any(t)
                X(:,t) = bsxfun(@times,X(:,t),(xi(:,k)==j));
            end
        end
    else
        % One row for each grid point; loop over powers of this predictor
        for j=1:max(psmatrix(k,:))
            t = (psmatrix(k,:)==j);
            X(:,t) = bsxfun(@times,X(:,t),xi(:,k).^j);
        end
    end
end

% Compute estimated fit and its variance
fxi = X*model.Coefs;
if nargout>=2
    fxiVar = X*model.CoefficientCovariance*X';
end


% ----------------------------------------------------------------------
function [xrow,psmatrix,psflag] = reduceterm(model,vnum,terminfo)
% Remove variables specified by vnum from all terms, and compute xrow as
% the mean of the remaining term. Also compute psmatrix as a matrix with
% one row per vnum value, specifying the corresponding variable's power
% (continuous) or setting (categorical), and psflag as a logical vector
% that is 1 if vnum is categorical.

xrow = zeros(size(terminfo.designTerms));
psmatrix = zeros(length(vnum),length(terminfo.designTerms));
psflag = terminfo.isCatVar(vnum);

for j=1:size(terminfo.terms,1)
    v = terminfo.terms(j,:);
    tj = terminfo.designTerms==j;
    pwr = v(vnum);
    meanx = gettermmean(v,vnum,model,terminfo);
    
    if all(pwr==0 | ~psflag)
        % Special case: removing term doesn't affect term size, because
        % vnum specifies only continous predictors and categorical ones
        % that are not part of this term
        xrow(tj) = meanx;
        psmatrix(:,tj) = repmat(pwr',1,sum(tj));
    elseif isscalar(vnum) && sum(terminfo.isCatVar(v>0))==sum(psflag)
        % Special case: vnum specifies the only categorical part of term
        xrow(tj) = meanx;
        psmatrix(:,tj) = 2:terminfo.numCatLevels(vnum);
    else
        % General case: suppose term is A*B*C*D*E and vnum is [D F B].
        % Here A represents both a categorical predictor and the number of
        % degrees of freedom for it.
        isreduced = ismember(find(v>0),vnum);
        termcatdims = terminfo.numCatLevels(v>0);           % size of ABCDE
        sz1 = ones(1,max(2,length(termcatdims)));
        sz1(~isreduced) = max(1,termcatdims(~isreduced)-1); % size of A1C1E
        sz2 = ones(1,max(2,length(termcatdims)));
        sz2(isreduced) = max(1,termcatdims(isreduced)-1);   % size of 1B1D1
        
        % Propagate mean of reduced term across the full term
        meanx = reshape(meanx,sz1);                  % shape to A1C1E
        meanx = repmat(meanx,sz2);                   % replicate to ABCDE
        xrow(tj) = meanx(:)';
        
        % Fill in powers for reduced continuous predictors
        controws = (pwr>0) & ~psflag;
        psmatrix(controws,tj) = repmat(pwr(controws),1,sum(tj));
        
        % Fill in settings for reduced categorical predictors (DB)
        catrows = (pwr>0) & psflag;
        catsettings = 1+fullfact(terminfo.numCatLevels(vnum(catrows))-1)';
        idx = reshape(1:size(catsettings,2),sz2);
        idx = repmat(idx,sz1);
        psmatrix(catrows,tj) = catsettings(:,idx(:));
    end
end


% ----------------------------------------------------------------------
function meanx = gettermmean(v,vnum,model,terminfo)
% Get the mean of the design matrix for a term after removing one or more
% variables from it

% Remove this variable from the term to get a subterm, for example remove B
% from A*B*C to get A*C
v(vnum) = 0;

% See if we have this term already, so we will already have means for it
[ok,row] = ismember(v,terminfo.terms,'rows');
if ok
    % Typically we have the subterm, so get the pre-computed means
    meanx = terminfo.designMeans(terminfo.designTerms==row);
else
    % We have to compute the design matrix columns for this term
    X = model.Data;
    if isstruct(X)
        X = X.X;      % use only the predictor data
        v(end) = [];  % omit response column from term
    end
    design = classreg.regr.modelutils.designmatrix(X,'Model',v,'VarNames',model.Formula.VariableNames);
    meanx = mean(design,1);
end
