
%   Copyright 2011 The MathWorks, Inc.


classdef QuadraticDiscriminantImpl < classreg.learning.impl.DiscriminantImpl
        
    properties(GetAccess=public,SetAccess=protected)
        % Predictor standard deviations per class, 1-by-p-by-K
        D = [];
        
        % Cell array with K diagonal S matrices from SVD for K classes.
        % Every S matrix is stored as a column-vector.
        S = {};
        
        % Cell array with K DiscriminantImpl objects. As of Aug 2011, they
        % are always LinearDiscriminantImpl objects.
        LinearPerClass = {};
    end
    
    properties(GetAccess=public,SetAccess=public,Hidden=true,Dependent=true)
        SaveMemory;
    end
    
    properties(GetAccess=public,SetAccess=protected,Hidden=true,Dependent=true)
        Sigma;
        InvSigma;
        LogDetSigma;
        MinGamma;
    end
        
    properties(Constant)
        AllowedTypes = {'quadratic' 'diagQuadratic' 'pseudoQuadratic'};
    end
    
    methods        
        function sm = get.SaveMemory(this)
            sm = this.LinearPerClass{1}.SaveMemory;
        end
        
        function this = set.SaveMemory(this,sm)
            K = numel(this.LinearPerClass);
            for k=1:K
                this.LinearPerClass{k}.SaveMemory = sm;
            end
        end
        
        function s = get.Sigma(this)
            K = numel(this.LinearPerClass);
            s1 = this.LinearPerClass{1}.Sigma;
            s = zeros([size(s1) K]);
            s(:,:,1) = s1;
            for k=2:K
                s(:,:,k) = this.LinearPerClass{k}.Sigma;
            end
        end
        
         function s = get.InvSigma(this)
            K = numel(this.LinearPerClass);
            s1 = this.LinearPerClass{1}.InvSigma;
            s = zeros([size(s1) K]);
            s(:,:,1) = s1;
            for k=2:K
                s(:,:,k) = this.LinearPerClass{k}.InvSigma;
            end
         end       
        
        function logdetsig = get.LogDetSigma(this)
            K = numel(this.LinearPerClass);
            logdetsig = zeros(K,1);
            for k=1:K
                logdetsig(k) = this.LinearPerClass{k}.LogDetSigma;
            end
        end
        
        function mingam = get.MinGamma(this)
            K = numel(this.LinearPerClass);
            mingam = zeros(K,1);
            for k=1:K
                mingam(k) = this.LinearPerClass{k}.MinGamma;
            end
            mingam = max(mingam);
        end
    end
        
    methods
        % X/sqrt(diag(Sigma))/Corr, N-by-p-by-K
        function m = linear(this,X)
            K = numel(this.LinearPerClass);
            m = zeros([size(X) K]);
            for k=1:K
                m(:,:,k) = linear(this.LinearPerClass{k},X);
            end
        end
        
        % X1*inv(Sigma)*X2', N-by-1-by-K
        function v = quadratic(this,X1,X2)
            K = numel(this.LinearPerClass);
            N = size(X1,1);
            v = zeros([N 1 K]);
            for k=1:K
                v(:,1,k) = quadratic(this.LinearPerClass{k},X1,X2);
            end
        end
        
        % Square of Mahalanobis distance for class indices in K,
        % N-by-numel(K)
        function m = mahal(this,K,X)
            m = zeros(size(X,1),numel(K));
            i = 0;
            for k=K
                i = i+1;
                m(:,i) = mahal(this.LinearPerClass{k},1,X);
            end
        end
        
        % Linear coefficients for classes i and j, p-by-1
        function v = linearCoeffs(this,i,j)
            mu1 = this.Mu(i,:);
            d1 = this.D(1,:,i);
            mu2 = this.Mu(j,:);
            d2 = this.D(1,:,j);
            v = bsxfun(@rdivide,linear(this.LinearPerClass{i},mu1),d1) ...
                - bsxfun(@rdivide,linear(this.LinearPerClass{j},mu2),d2);
            v = v';
        end
        
        % Quadratic coefficients for classes i and j, p-by-p
        function m = quadraticCoeffs(this,i,j)
            m = -0.5*( this.LinearPerClass{i}.InvSigma - this.LinearPerClass{j}.InvSigma );
        end
        
        % Const term. For QDA we must take the difference of LogDetSigma's
        % because the covariance matrices are not equal.
        function c = constantTerm(this,i,j)
            p = size(this.Mu,2);
            c = 0.5*( mahal(this.LinearPerClass{i},1,zeros(1,p)) ...
                - mahal(this.LinearPerClass{j},1,zeros(1,p)) );
            c = c + 0.5*( this.LinearPerClass{i}.LogDetSigma - this.LinearPerClass{j}.LogDetSigma );
        end

        function delran = deltaRange(this,gamma)
            delran = zeros(1,2);
        end
        
        function delpred = deltaPredictor(this,gamma)
            delpred = zeros(1,size(this.Mu,2));
        end
    
        function nCoeffs = nLinearCoeffs(this,delta)
            if ~isnumeric(delta) || ~isvector(delta)
                error(message('stats:classreg:learning:impl:QuadraticDiscriminantImpl:nLinearCoeffs:BadDelta'));
            end
            nCoeffs = repmat(size(this.Mu,2),1,numel(delta));
        end
        
        function this = QuadraticDiscriminantImpl(type,d,s,v,gamma,delta,mu,classWeights,saveMemory)
            % Disallow delta
            if delta~=0
                error(message('stats:classreg:learning:impl:QuadraticDiscriminantImpl:QuadraticDiscriminantImpl:DeltaNotAllowed'));
            end
            
            this = this@classreg.learning.impl.DiscriminantImpl(gamma,0,mu,classWeights);            
            this.D = d;
            this.S = s;
            
            type = spellCompleteType(type);

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

            type = convertTypeToLinear(type);
            
            K = size(mu,1);
            this.LinearPerClass = cell(K,1);
            for k=1:K
                this.LinearPerClass{k} = classreg.learning.impl.LinearDiscriminantImpl(...
                    type,d(1,:,k),s{k},v{k},gamma,0,mu(k,:),classWeights,saveMemory);
            end
        end
                
        function this = setType(this,type)
            type = spellCompleteType(type);
            if ~strcmpi(type,this.Type)
                d = this.D;
                s = this.S;
                [~,~,gamma] = forceType(type,d,s);
                this.Type = type;
                this.Gamma = gamma;
                
                type = convertTypeToLinear(type);
                K = numel(this.LinearPerClass);
                for k=1:K
                    this.LinearPerClass{k} = setTypeAndGamma(this.LinearPerClass{k},type,gamma);
                end
            end
        end
        
        function this = setGamma(this,gamma)
            if gamma~=this.Gamma
                d = this.D;
                s = this.S;
                allowTypes = allowedTypes(d,s,gamma);
                if ~ismember(this.Type,allowTypes)
                    this.Type = forceGamma(size(this.Mu,2),s,gamma);
                end
                
                type = convertTypeToLinear(this.Type);
                K = numel(this.LinearPerClass);
                for k=1:K
                    this.LinearPerClass{k} = setTypeAndGamma(this.LinearPerClass{k},type,gamma);
                end
                this.Gamma = gamma;
            end
        end
        
        function this = setDelta(this,delta)
            if delta>0
                error(message('stats:classreg:learning:impl:QuadraticDiscriminantImpl:setDelta:Disallowed'));
            end
        end
    end
    
end


function mingam = minGamma(p,s)
K = numel(s);
mingam = zeros(K,1);
for k=1:K
    mingam(k) = classreg.learning.impl.LinearDiscriminantImpl.minGamma(p,s{k});
end
mingam = max(mingam);
end

function tf = canBeFull(d,s,gamma)
K = numel(s);
tf = false(K,1);
for k=1:K
    tf(k) = classreg.learning.impl.LinearDiscriminantImpl.canBeFull(d(1,:,k),s{k},gamma);
end
tf = all(tf);
end

function tf = canBePseudo(d,s,gamma)
tf = gamma==0;
end

function tf = canBeDiagonal(d,s,gamma)
tf = gamma==1;
end

function [d,s,gamma] = forceFull(d,s)
gamma = minGamma(size(d,2),s);
if gamma>0
    error(message('stats:classreg:learning:impl:QuadraticDiscriminantImpl:forceFull:ZeroGammaDisallowed'));
end
end

function [d,s,gamma] = forcePseudo(d,s)
gamma = 0;
end

function [d,s,gamma] = forceDiagonal(d,s)
gamma = 1;
end

function [d,s,gamma] = forceType(type,d,s)
switch lower(type)
    case 'quadratic'
        [d,s,gamma] = forceFull(d,s);
    case 'pseudoquadratic'
        [d,s,gamma] = forcePseudo(d,s);
    case 'diagquadratic'
        [d,s,gamma] = forceDiagonal(d,s);
end
end

function type = forceGamma(p,s,gamma)
mingam = minGamma(p,s);
if     gamma==1
    type = 'diagQuadratic';
elseif gamma<mingam
    error(message('stats:classreg:learning:impl:QuadraticDiscriminantImpl:forceGamma:GammaTooSmall', sprintf( '%g', mingam )));
else
    type = 'quadratic';
end
end

function types = allowedTypes(d,s,gamma)
if gamma>0 && gamma<1
    error(message('stats:classreg:learning:impl:QuadraticDiscriminantImpl:allowedTypes:GammaBetween0and1Disallowed'));
end
types = {};
if     canBeFull(d,s,gamma)
    types = {'quadratic'};
end
if canBePseudo(d,s,gamma)
    types = [types {'pseudoQuadratic'}];
end
if canBeDiagonal(d,s,gamma)
    types = [types {'diagQuadratic'}];
end
end

function type = spellCompleteType(type)
allowedTypes = classreg.learning.impl.QuadraticDiscriminantImpl.AllowedTypes;
tf = strncmpi(type,allowedTypes,length(type));
if sum(tf)~=1
    error(message('stats:classreg:learning:impl:QuadraticDiscriminantImpl:spellCompleteType:BadType', sprintf( ' %s', allowedTypes{ : } )));
end
type = allowedTypes{tf};
end

function type = convertTypeToLinear(type)
if ~isempty(type)
    switch lower(type)
        case 'diagquadratic'
            type = 'diagLinear';
        case 'pseudoquadratic'
            type = 'pseudoLinear';
        case 'quadratic'
            type = 'linear';
        otherwise
            error(message('stats:classreg:learning:impl:QuadraticDiscriminantImpl:convertTypeToLinear:BadType', type));
    end
end
end
