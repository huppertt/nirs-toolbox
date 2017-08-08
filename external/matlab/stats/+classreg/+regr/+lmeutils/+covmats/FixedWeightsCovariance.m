classdef FixedWeightsCovariance < classreg.regr.lmeutils.covmats.CovarianceMatrix
%  This feature is intended for internal use only and is subject to change
%  at any time without warning. 

%FixedWeightsCovariance - Class to represent a 'Fixed Weights' covariance matrix. 
%   This is a helper class for LinearMixedModel. Suppose D(theta) is a 
%   matrix that is parameterized by an unconstrained vector theta. The
%   covariance matrix PSI represented by this class can be written like
%   this:
% 
%          PSI = sigma^2 * D(theta) = sigma^2 * L(theta)*L(theta)'
%   where
%          PSI = covariance matrix (positive definite and symmetric)
%            D = scaled covariance matrix (positive definite and symmetric)
%            L = lower triangular Cholesky factor of D
%        theta = unconstrained parameter vector for D
%      sigma^2 = positive multiplier of D in the definition of PSI
%
%   When PSI is of Type 'FixedWeights', only the diagonal elements of PSI 
%   are non-zero, all other elements are zero. In this case D is a diagonal
%   matrix with known diagonal elements. Here's an example to illustrate 
%   what the unconstrained, natural and canonical parameterizations look 
%   like in this case. Suppose W is a q by 1 vector and suppose elements 
%   1/W(i) are arranged along the diagonal of the q by q matrix D:
%
%   Suppose D = [1/W(1)   0      0
%                 0     1/W(2)   0
%                ..      ..     ..
%                 0       0    1/W(q)]
%
%   The positive definiteness of D requires that 1/W(i) > 0 which implies
%   that W(i) is any finite positive value (including 0). All elements of D
%   are known and there are no free parameters. Therefore the unconstrained
%   parameters are given by:
%  
%           theta = [];
%
%   There are also no Natural parameters and hence:
%
%           eta = [];
%
%   There are no Canonical parameters:
%
%           heta = [];
%
%   FixedWeightsCovariance methods:
%       FixedWeightsCovariance - Create a 'Fixed Weights' covariance matrix.
%       getUnconstrainedParameters - Get unconstrained parameters
%       getNaturalParameters - Get natural parameters
%       getCanonicalParameters - Get canonical parameters
%       setUnconstrainedParameters - Set unconstrained parameters
%       setNaturalParameters - Set natural parameters
%       getLowerTriangularCholeskyFactor - Get lower triangular Cholesky factor
%       getSigma - Get sigma
%       createUnconstrainedParameterization - Create unconstrained parameterization
%       createNaturalParameterization - Create Natural parameterization
%       createCanonicalParameterization - Create Canonical parameterization
%
%   FixedWeightsCovariance properties:
%       Name - Name of this covariance matrix
%       Size - Size of this covariance matrix
%       VariableNames - Names of random variables whose covariance is PSI
%       Type - The type of this covariance matrix (e.g., "Diagonal").
%       CovariancePattern - Pattern of non-zeros in this covariance matrix
%       WhichParameterization - Parameterization that is currently in use.
%       NumParametersExcludingSigma - Number of parameters excluding sigma.
%       Weights - A vector of fixed weights used to create this covariance.
%
%   See also classreg.regr.lmeutils.covmats.CovarianceMatrix
    
%   Copyright 2012-2014 The MathWorks, Inc. 

    properties (Access=protected)        
%An indicator matrix for the diagonal elements. Basically this is equal to
%   logical(eye(this.Size)).        
        DiagonalElements
        
%Weights - A vector of length Size by 1 used to create this covariance.
        Weights        
    end
   
    methods (Access=public)
        function this = FixedWeightsCovariance(dimension,varargin)
%FixedWeightsCovariance - Construct a FixedWeightsCovariance object.
%   OBJ = FixedWeightsCovariance(DIMENSION) takes an integer DIMENSION and 
%   constructs a FixedWeightsCovariance object.
%
%   OBJ = FixedWeightsCovariance(DIMENSION,'VariableNames',VNAMES) also 
%   specifies the variable names associated with rows and columns of this 
%   covariance matrix.
%
%   OBJ = FixedWeightsCovariance(DIMENSION,..,'Name',NAME) also specifies 
%   the name to be associated with this covariance matrix.
%
%   OBJ = FixedWeightsCovariance(DIMENSION,...,'Weights',WEIGHTS) also 
%   specifies a vector WEIGHTS of length DIMENSION that is to be used in 
%   constructing the FixedWeightsCovariance. 
        
            this = this@classreg.regr.lmeutils.covmats.CovarianceMatrix(dimension,varargin{:});
            this.DiagonalElements = logical(eye(this.Size));

            % Parse weights.
            dfltWeights = [];
            parnames = {'Weights'};
            dflts = {dfltWeights};
            [w,~,~] = internal.stats.parseArgs(parnames,dflts,varargin{:});

            if isempty(w)
                w = ones(dimension,1);
            else
                % w = vector, real, >= 0 and not Inf.
                if ~isvector(w) || ~isnumeric(w) || ~isreal(w) || ~all(w >= 0) || any(isinf(w))
                    % <entry key="BadWeights">WEIGHTS must be a finite, real and non-negative vector.</entry>
                    error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadWeights'));
                end
                
                % Make w into a column vector.
                if size(w,1) == 1
                    w = w';
                end
                
                if length(w) ~= dimension
                    % <entry key="BadWeightsLength">WEIGHTS must be a vector of length {0}.</entry>
                    error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadWeightsLength',num2str(dimension)));
                end                
            end

            % Set weights.
            this.Weights = w;
        
        end % end of FixedWeightsCovariance. 
        
        function this = setWeights(this,w)
%setWeights - Set the weight vector for this FixedWeightsCovariance.
%   OBJ = setWeights(OBJ,NEWWEIGHTS) takes a finite, real and non-negative
%   weight vector NEWWEIGHTS and returns a FixedWeightsCovariance object
%   with the new weights.

            % w = vector, real, >= 0 and not Inf.
            if ~isvector(w) || ~isnumeric(w) || ~isreal(w) || ~all(w >= 0) || any(isinf(w))
                % <entry key="BadWeights">WEIGHTS must be a finite, real and non-negative vector.</entry>
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadWeights'));
            end

            % Make w into a column vector.
            if size(w,1) == 1
                w = w';
            end

            if length(w) ~= this.Size
                % <entry key="BadWeightsLength">WEIGHTS must be a vector of length {0}.</entry>
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadWeightsLength',num2str(this.Size)));
            end                

            % Set weights.
            this.Weights = w;
            
        end % end of setWeights.
        
        function w = getWeights(this)
%getWeights - Get the weight vector for this FixedWeightsCovariance.
%   WEIGHTS = getWeights(OBJ) returns the weight vector for the object OBJ
%   representing a FixedWeightsCovariance.
            
            w = this.Weights;

        end % end of getWeights.        
        
        function verifyTransformations(this)
%verifyTransformations - Check the parameter transformations between 
%   unconstrained, natural and canonical for accuracy.
          
            % (1) theta to D and D to theta
            theta = randn(0,1);
            D = theta2D(this,theta);            
            theta_recon = D2theta(this,D);            
            if isempty( max(abs(theta(:) - theta_recon(:))) )
                err(1) = 0;
            else
                err(1) = 1;
            end

            % (2) D to theta and theta to D
            weights = rand(this.Size,1);
            this = setWeights(this,weights);
            D = diag(1./weights);
            L = this.getLowerTriangularCholeskyFactor;
            D_recon = L*L';           
            err(2) = max(abs(D(:) - D_recon(:)));
                                    
            % (3) 
            % (A) L via chol vs. getLowerTriangularCholeskyFactor
            theta = randn(0,1);
            D = theta2D(this,theta);            
            L = diag(sqrt(diag(D)));            
            this = setUnconstrainedParameters(this,theta);
            L1 = getLowerTriangularCholeskyFactor(this);            
            err(3) = max(abs(L(:) - L1(:)));

            % (4)
            % (B) theta to natural and then getLowerTriangularCholeskyFactor
            natural = convertUnconstrainedToNatural(this,theta);
            this = setNaturalParameters(this,natural);
            L2 = getLowerTriangularCholeskyFactor(this);            
            err(4) = max(abs(L(:) - L2(:)));
            
            % (5) 
            % (C) canonical parameters from unconstrained vs. natural
            this = setUnconstrainedParameters(this,theta);
            canonical1 = getCanonicalParameters(this);
            this = setNaturalParameters(this,natural);
            canonical2 = getCanonicalParameters(this);
            if isempty( max(abs(canonical1(:) - canonical2(:))) )
                err(5) = 0;
            else
                err(5) = 1;
            end
                                    
            if max(err(:)) <= sqrt(eps)
                disp(getString(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:String_TransformsOK')));
            else
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:String_TransformsNotOK'));
            end            
            
        end % end of verifyTransformations.
                
    end % end of verifyTransformations.
    
     methods (Access=public, Hidden=true)
      
        % A subclass must know how to initialize theta and how to convert
        % to L. For verifying transforms, it is also useful for subclasses
        % to define theta2D and D2theta methods.
        function theta = initializeTheta(this) %#ok<MANU>
        % Return empty theta.   
            
            theta = zeros(0,1);
            
        end % end of initializeTheta.
        
        function L = theta2L(this,theta)
        % L(i,i) = 1/sqrt(W(i)). Return sparse L.    
            
            assert( isempty(theta) );
            q = this.Size;
            L = spdiags(1./sqrt(this.Weights),0,q,q);
            
        end % end of theta2L.
        
        function D = theta2D(this,theta)
        % D(i,i) = 1/W(i). Return sparse D.  
        
            assert( isempty(theta) );
            q = this.Size;
            D = spdiags(1./this.Weights,0,q,q);
            
        end % end of theta2D.
        
        function theta = D2theta(this,D) %#ok<INUSD>
        % For a FixedWeightsCovariance, theta = [];    
            
            theta = zeros(0,1);
        
        end % end of D2theta.
               
        % A subclass must know how to convert between unconstrained and
        % natural parameterizations.
        function natural = convertUnconstrainedToNatural(this,unconstrained) %#ok<INUSL>
        % We assert that unconstrained is [] and then return natural = [].
        
            assert( isempty(unconstrained) );
            natural = zeros(0,1);
            
        end % end of convertUnconstrainedToNatural.
        
        function unconstrained = convertNaturalToUnconstrained(this,natural) %#ok<INUSL>
        % We assert that natural is [] and then return unconstrained = [].
        
            assert( isempty(natural) );            
            unconstrained = zeros(0,1);
            
        end % end of convertNaturalToUnconstrained.
        
        % A subclass must know how to convert between natural and canonical
        % parameterizations.
        function canonical = convertNaturalToCanonical(this,natural) %#ok<INUSL>
        % We assert that natural is [] and then return canonical = [].
        
            assert( isempty(natural) );            
            canonical = zeros(0,1);
            
        end % end of convertNaturalToCanonical.
        
        function natural = convertCanonicalToNatural(this,canonical) %#ok<INUSL>
        % We assert that canonical is [] and then return natural = [].    
            
            assert( isempty(canonical) );            
            natural = zeros(0,1);
        
        end % end of convertCanonicalToNatural.
        
        % A subclass must define canonical parameter names.
        function names = defineCanonicalParameterNames(this) %#ok<MANU>
        % Since there are no canonical parameters, names is an empty table.
            names = table({},'VariableNames',{'Name'});
        end % end of defineCanonicalParameterNames.
        
        % A subclass must define its type (e.g., 'Full','Diagonal' etc.)
        function type = defineType(this) %#ok<MANU>
        % A 'FixedWeightsCovariance' has type 'Fixed Weights'.
            type = classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FIXEDWEIGHTS;        
        end % end of defineType.
        
        % A subclass must define its covariance pattern.
        function covariancepattern = defineCovariancePattern(this)
        % A FixedWeightsCovariance has non-zeros only on the diagonal 
        % elements.   
            covariancepattern = logical(eye(this.Size));
            
        end % end of defineCovariancePattern.
        
    end
        
end