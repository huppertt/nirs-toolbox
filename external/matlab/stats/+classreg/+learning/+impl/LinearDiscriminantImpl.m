
%   Copyright 2011 The MathWorks, Inc.


classdef LinearDiscriminantImpl < classreg.learning.impl.DiscriminantImpl
    
    % Factorization for the covariance matrix:
    % D = std(X); [~,S,V] = svd(X*diag(1./D),'econ')
    %
    % For non-singular Sigma:
    % R = diag(S)*V'*diag(D); Sigma = R'*R
    %
    % By convention, D, S and V are never changed by the class. If you need
    % to adjust them to account for various conditions on Gamma and Type
    % before passing into makeCalculator, you need to make local copies.
    properties(GetAccess=public,SetAccess=protected)
        D = [];
        S = [];
        V = [];
    end
    
    properties(GetAccess=protected,SetAccess=protected)
        PrivCalculator = [];
    end
    
    properties(GetAccess=public,SetAccess=public,Hidden=true,Dependent=true)
        SaveMemory;
    end

    properties(GetAccess=public,SetAccess=protected,Hidden=true,Dependent=true)
        Calculator;
        Sigma;
        InvSigma;
        LogDetSigma;
        MinGamma;
    end
        
    properties(Constant)
        AllowedTypes = {'linear' 'diagLinear' 'pseudoLinear'};
    end
    
    methods
        function sm = get.SaveMemory(this)
            sm = isempty(this.PrivCalculator);
        end
        
        function this = set.SaveMemory(this,sm)
            if sm && ~isempty(this.PrivCalculator)
                this.PrivCalculator = [];
            end
            if ~sm && isempty(this.PrivCalculator)
                this.PrivCalculator = fetchCalculator(this);
            end
        end
        
        function calc = get.Calculator(this)
            if isempty(this.PrivCalculator)
                calc = fetchCalculator(this);
            else
                calc = this.PrivCalculator;
            end
        end
        
        function s = get.Sigma(this)
            s = sigma(this.Calculator,this.D,this.S,this.V);
        end
        
        function s = get.InvSigma(this)
            s = invSigma(this.Calculator);
        end
        
        function logdetsig = get.LogDetSigma(this)
            logdetsig = logDetSigma(this.Calculator,this.D,this.S,this.V);
        end
        
        function mg = get.MinGamma(this)
            mg = classreg.learning.impl.LinearDiscriminantImpl.minGamma(size(this.Mu,2),this.S);
        end
    end
        
    methods
        % X/sqrt(diag(Sigma))/Corr, N-by-p
        function m = linear(this,X)
            m = linear(this.Calculator,X);
        end
        
        % X1*inv(Sigma)*X2', N-by-1        
        function v = quadratic(this,X1,X2)
            v = quadratic(this.Calculator,X1,X2);
        end
        
        % Square of Mahalanobis distance for class indices in K,
        % N-by-numel(K)
        function m = mahal(this,K,X)
            m = mahal(this.Calculator,K,X);
        end
        
        % Linear coefficients for classes i and j, p-by-1
        function v = linearCoeffs(this,i,j)
            v = linearCoeffs(this.Calculator,i,j);
        end
        
        % Const term
        function c = constantTerm(this,i,j)
            c = constantTerm(this.Calculator,i,j);
        end

        function delran = deltaRange(this,gamma)
            if nargin<2
                gamma = this.Gamma;
            end
            r = deltaPredictor(this,gamma);
            delran = [min(r(:)) max(r(:))];
        end
        
        function delpred = deltaPredictor(this,gamma)
            if nargin<2 || gamma==this.Gamma
                calc = this.Calculator;
                if isa(calc,'classreg.learning.impl.RegularizedDiscriminantCalculator')
                    delpred = max(abs(calc.CenteredScaledMuOverCorr));
                else
                    delpred = max(abs(linear(calc,this.CenteredMu)));
                end
            else
                calc = makeCalculator(...
                    this.D,this.S,this.V,gamma,0,this.Mu,this.BetweenMu,this.CenteredMu);
                delpred = max(abs(linear(calc,this.CenteredMu)));
            end
        end
        
        function nCoeffs = nLinearCoeffs(this,delta)
            if ~isnumeric(delta) || ~isvector(delta)
                error(message('stats:classreg:learning:impl:LinearDiscriminantImpl:nLinearCoeffs:BadDelta'));
            end
            delta = delta(:);
            
            Ndelta = numel(delta);
            
            p = size(this.Mu,2);
            if any(delta>0)
                maxDeltaPerPredictor = deltaPredictor(this);
                nCoeffs = sum(bsxfun(@ge,maxDeltaPerPredictor,repmat(delta,1,p)),2);
            else
                nCoeffs = repmat(p,1,Ndelta);
            end
        end

        function this = LinearDiscriminantImpl(type,d,s,v,gamma,delta,mu,classWeights,saveMemory)
            this = this@classreg.learning.impl.DiscriminantImpl(gamma,delta,mu,classWeights);
            this.D = d;
            this.S = s;
            this.V = v;
            
            type = spellCompleteType(type);
            
            % If gamma==0, consider this as a sign that the user does not
            % know how to regularize the data. In this case, set gamma to
            % the correct value based on the type of the discriminant.
            % Otherwise make sure the value of gamma is consistent with the
            % discriminant type.
            if gamma==0
                [d,s,gamma] = forceType(type,d,s);
            else
                allowTypes = allowedTypes(d,s,gamma);
                if ~ismember(lower(type),lower(allowTypes))
                    type = forceGamma(size(mu,2),s,gamma);
                end
            end
            
            this.Type = type;
            this.Gamma = gamma;

            % Make a calculator object unless SaveMemory is on. The
            % calculator object stores the full covariance matrix.
            if ~saveMemory
                this.PrivCalculator = makeCalculator(...
                    d,s,v,gamma,delta,mu,this.BetweenMu,this.CenteredMu);
            end
        end
        
        function this = setType(this,type)
            type = spellCompleteType(type);
            
            % Reset type only if it is different from the one currently in
            % use.
            if ~strcmpi(type,this.Type)
                d = this.D;
                s = this.S;
                
                % Force gamma to comply with the new type
                [d,s,gamma] = forceType(type,d,s);
                this.Type = type;

                % If the calculator object exists (meaning SaveMemory is
                % off), make a new one for the new value of gamma
                if this.Gamma~=gamma
                    this.Gamma = gamma;                    
                    if ~isempty(this.PrivCalculator)
                        this.PrivCalculator = makeCalculator(...
                            d,s,this.V,gamma,this.Delta,this.Mu,this.BetweenMu,this.CenteredMu);
                    end
                end
            end
        end
                
        function this = setGamma(this,gamma)
            % Reset gamma only if it different from the one in use. Beware
            % of an edge case: type=pseudoLinear for a singular covariance
            % matrix. In this case, gamma is forced to be 0. If the user
            % assigns 0, no action is taken because of the check below. But
            % if the user assigns gamma between 0 and MinGamma, he gets an
            % error message.
            if gamma~=this.Gamma
                d = this.D;
                s = this.S;
                
                % What discriminant types are allowed for this gamma? If
                % the current type is not one of them, force it to be the
                % right type for this gamma
                allowTypes = allowedTypes(d,s,gamma);
                if ~ismember(this.Type,allowTypes)
                    this.Type = forceGamma(size(this.Mu,2),s,gamma);
                end
                
                this.Gamma = gamma;
                
                % If the calculator object exists (meaning SaveMemory is
                % off), make a new one for the new value of gamma
                if ~isempty(this.PrivCalculator)
                    [d,s] = forceType(this.Type,d,s);                
                    this.PrivCalculator = makeCalculator(...
                        d,s,this.V,gamma,this.Delta,this.Mu,this.BetweenMu,this.CenteredMu);
                end
            end
        end
        
        % Set type and gamma simultaneously. As of Aug 2011, this method is
        % only called by QuadraticDiscriminantImpl which first finds type
        % and gamma appropriate for all classes and then loops over all
        % linear discriminants resetting them in turn.
        %
        % This method is used to avoid constructing a new calculator object
        % both in setType and setGamma methods because this is potentially
        % an expensive operation. First we empty PrivCalculator which
        % ensures that the new calculator is not constructed. Then we
        % construct it once, if needed, after gamma and delta are set to
        % the new values.
        function this = setTypeAndGamma(this,type,gamma)
            oldCalculator = this.PrivCalculator;
            this.PrivCalculator = [];
            this = setType(this,type);
            this = setGamma(this,gamma);
            if ~isempty(oldCalculator)
                d = this.D;
                s = this.S;
                [d,s] = forceType(this.Type,d,s);
                this.PrivCalculator = makeCalculator(...
                        d,s,this.V,this.Gamma,this.Delta,this.Mu,this.BetweenMu,this.CenteredMu);
            end
        end
        
        function this = setDelta(this,delta)
            if delta~=this.Delta
                this.Delta = delta;
                if ~isempty(this.PrivCalculator)
                    % If the calculator is already
                    % RegularizedDiscriminantCalculator, no need to
                    % recompute anything; just assign the new delta to its
                    % property. Otherwise make a new calculator.
                    if delta>0 && isa(this.PrivCalculator,...
                            'classreg.learning.impl.RegularizedDiscriminantCalculator')
                        this.PrivCalculator.Delta = delta;
                    else
                        this.PrivCalculator = ...
                            makeCalculator(...
                            this.D,this.S,this.V,this.Gamma,delta,...
                            this.Mu,this.BetweenMu,this.CenteredMu);
                    end
                end                
            end
        end
    end
    
    methods(Access=protected)
        % Implements get.Calculator method. The covariance matrix is
        % computed on the fly every time it is asked for.
        function calc = fetchCalculator(this)
            d = this.D;
            s = this.S;
            v = this.V;

            [d,s] = forceType(this.Type,d,s);
            gamma = this.Gamma;
            
            calc = makeCalculator(...
                d,s,v,gamma,this.Delta,this.Mu,this.BetweenMu,this.CenteredMu);
        end
    end
    
    methods(Static,Hidden)
        function mingam = minGamma(p,s)
            s = s.^2;
            maxs = max(s);
            if any(s < p*eps(maxs))
                mingam = p*eps(maxs);
            else
                mingam = 0;
            end
        end
        
        function tf = canBeFull(d,s,gamma)
            tf = all(d>0) && gamma<1 && ...
                (all(s>0) || gamma>classreg.learning.impl.LinearDiscriminantImpl.minGamma(numel(d),s));
        end
        
        function tf = canBePseudo(d,s,gamma)
            tf = gamma==0;
        end
        
        function tf = canBeDiagonal(d,s,gamma)
            tf = gamma==1;
        end
    end
    
end


function [d,s,gamma] = forceFull(d,s)
gamma = classreg.learning.impl.LinearDiscriminantImpl.minGamma(numel(d),s);
end

function [d,s,gamma] = forcePseudo(d,s)
% Set zero elements in D and S to Inf, so their inverse values
% would be 0.
d(d==0) = Inf;
s(s==0) = Inf;
gamma = 0;
end

function [d,s,gamma] = forceDiagonal(d,s)
% Set zero elements in D to Inf, so its inverse values would be
% 0.
d(d==0) = Inf;
gamma = 1;
end

function [d,s,gamma] = forceType(type,d,s)
switch lower(type)
    case 'linear'
        [d,s,gamma] = forceFull(d,s);
    case 'pseudolinear'
        [d,s,gamma] = forcePseudo(d,s);
    case 'diaglinear'
        [d,s,gamma] = forceDiagonal(d,s);
end
end

function type = forceGamma(p,s,gamma)
mingam = classreg.learning.impl.LinearDiscriminantImpl.minGamma(p,s);
if     gamma==1
    type = 'diagLinear';
elseif gamma<mingam
    error(message('stats:classreg:learning:impl:LinearDiscriminantImpl:forceGamma:GammaTooSmall', sprintf( '%g', mingam )));
else
    type = 'linear';
end
end

function types = allowedTypes(d,s,gamma)
types = {};
if     classreg.learning.impl.LinearDiscriminantImpl.canBeFull(d,s,gamma)
    types = {'linear'};
end
if classreg.learning.impl.LinearDiscriminantImpl.canBePseudo(d,s,gamma)
    types = [types {'pseudoLinear'}];
end
if classreg.learning.impl.LinearDiscriminantImpl.canBeDiagonal(d,s,gamma)
    types = [types {'diagLinear'}];
end
end

function type = spellCompleteType(type)
allowedTypes = classreg.learning.impl.LinearDiscriminantImpl.AllowedTypes;
tf = strncmpi(type,allowedTypes,length(type));
if sum(tf)~=1
    error(message('stats:classreg:learning:impl:LinearDiscriminantImpl:spellCompleteType:BadType', sprintf( ' %s', allowedTypes{ : } )));
end
type = allowedTypes{tf};
end

% Make a calculator object based on values of gamma and delta. If any of
% the two regularization parameters is above 0, makeCalculator returns
% RegularizedDiscriminantCalculator. One and only exception is the case
% gamma=1 and delta=0 when DiagonalDiscriminantCalculator is made. The
% non-regularized CholeskyDiscriminantCalculator is returned if gamma=0 and
% delta=0.
function calc = makeCalculator(d,s,v,gamma,delta,mu,betweenMu,centeredMu)
if     gamma==1
    if delta>0
        calc = classreg.learning.impl.RegularizedDiscriminantCalculator(...
            mu,1./d,eye(size(mu,2)),gamma,delta,betweenMu,centeredMu);
    else
        calc = classreg.learning.impl.DiagonalDiscriminantCalculator(mu,1./d);
    end
elseif gamma>0
    p = size(mu,2);
    s2 = (1-gamma)*s.^2 + gamma;
    if any(abs(s2)<p*eps(max(abs(s2))))
        error(message('stats:classreg:learning:impl:LinearDiscriminantImpl:makeCalculator:PooledSingular'));
    end
    invCorr = bsxfun(@times,v,1./s2'-1/gamma)*v';
    invCorr(1:p+1:end) = invCorr(1:p+1:end) + 1/gamma;
    calc = classreg.learning.impl.RegularizedDiscriminantCalculator(...
        mu,1./d,invCorr,gamma,delta,betweenMu,centeredMu);
elseif delta>0 % gamma=0
    invR = bsxfun(@times,v,1./s');
    invCorr = invR*invR';
    calc = classreg.learning.impl.RegularizedDiscriminantCalculator(...
        mu,1./d,invCorr,gamma,delta,betweenMu,centeredMu);
else
    invR = bsxfun(@times,v,1./s');
    calc = classreg.learning.impl.CholeskyDiscriminantCalculator(...
        mu,1./d,invR);
end
end
