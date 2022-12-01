function Objective = statmvnrobj(Data, Design, Param, Covar, Resid, CholCovar, isdiagonal)
%STATMVNROBJ Log-likelihood for multivariate normal regression without missing data.

%    Copyright 2006-2007 The MathWorks, Inc.


Param = Param(:);

[NumSamples, NumSeries] = size(Data);

if nargin<6 || isempty(CholCovar)
    if (nargin>=7) && isdiagonal    % use faster formula for diagonal covariance
        d = diag(Covar);
        CholState = sum(d<=0);
        CholCovar = diag(sqrt(d));
    else
        [CholCovar, CholState] = chol(Covar);
    end
else
    CholState = size(Covar,1) - rank(CholCovar);
end
if CholState > 0
    error(message('stats:mvregresslike:NonPosDefCov'));
end

LogTwoPi = log(2.0 * pi);
LogDetCovar = 2.0 * sum(log(diag(CholCovar)));

Count = 0;
Objective = 0.0;

if nargin>=5 && ~isempty(Resid)
    Objective = -0.5 * sum(sum((Resid / CholCovar).^2,2));
    Count = size(Resid,1);
elseif iscell(Design)
    if numel(Design) > 1
        for k = 1:NumSamples
            if ~any(isnan(Data(k,:)))
                Count = Count + 1;
                Resid = CholCovar' \ (Data(k,:)' - Design{k} * Param);
                Objective = Objective - 0.5 * (Resid' * Resid);
            end
        end
    else
        for k = 1:NumSamples
            if ~any(isnan(Data(k,:)))
                Count = Count + 1;
                Resid = CholCovar' \ (Data(k,:)' - Design{1} * Param);
                Objective = Objective - 0.5 * (Resid' * Resid);
            end
        end
    end
else
    for k = 1:NumSamples
        if ~isnan(Data(k))
            Count = Count + 1;
            Resid = CholCovar' \ (Data(k) - Design(k,:) * Param);
            Objective = Objective - 0.5 * (Resid' * Resid);
        end
    end
end

Objective = Objective - 0.5 * Count * (NumSeries * LogTwoPi + LogDetCovar);

if Count < 1
    Objective = NaN;
    warning(message('stats:mvregresslike:AllNaNData'));
end
