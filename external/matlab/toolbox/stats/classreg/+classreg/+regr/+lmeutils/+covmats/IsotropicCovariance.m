classdef IsotropicCovariance < classreg.regr.lmeutils.covmats.CovarianceMatrix
%  This feature is intended for internal use only and is subject to change
%  at any time without warning. 

%IsotropicCovariance - Class to represent a 'IsotropicCovariance' covariance matrix. 
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
%   When PSI is of Type 'Isotropic', only the diagonal elements of PSI are
%   non-zero, all other elements are zero and in addition all diagonal 
%   elements are equal. Here's an example to illustrate what the 
%   unconstrained, natural and canonical parameterizations look like in 
%   this case:
%
%   Suppose D = [D_11   0     0
%                 0    D_11   0
%                 0     0    D_11]
%
%   The positive definiteness of D requires that D_11 >= 0. An
%   unconstrained parameterization of D that guarantees positive
%   definiteness is given by:
%
%           theta = [log(D_11)] 
%
%   Suppose StdCorr is a matrix such that:
%              
%                StdCorr(i,i) = sqrt(PSI(i,i))
%                StdCorr(i,j) = 0
%
%   and suppose Omega is a matrix such that:
%
%               Omega(i,i) = log( StdCorr(i,i) )
%               Omega(i,j) = 0
%
%   then the Natural parameterization eta is given by the first row and 
%   first column element of Omega and the Canonical parameterization heta 
%   is given by the first row and first column element of StdCorr.
%
%   IsotropicCovariance methods:
%       IsotropicCovariance - Create a 'Isotropic' covariance matrix.
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
%   IsotropicCovariance properties:
%       Name - Name of this covariance matrix
%       Size - Size of this covariance matrix
%       VariableNames - Names of random variables whose covariance is PSI
%       Type - The type of this covariance matrix (e.g., "Diagonal").
%       CovariancePattern - Pattern of non-zeros in this covariance matrix
%       WhichParameterization - Parameterization that is currently in use.
%       NumParametersExcludingSigma - Number of parameters excluding sigma.    
%
%   See also classreg.regr.lmeutils.covmats.CovarianceMatrix
    
%   Copyright 2012-2013 The MathWorks, Inc.

    properties (Access=protected)        
%An indicator matrix for the diagonal elements. Basically this is equal to
%   logical(eye(this.Size)).        
        DiagonalElements
        
    end

    methods (Access=public)
        
         function this = IsotropicCovariance(dimension,varargin)
%IsotropicCovariance - Construct a IsotropicCovariance object.
%   OBJ = IsotropicCovariance(DIMENSION) takes an integer DIMENSION and 
%   constructs a IsotropicCovariance object.
%
%   OBJ = IsotropicCovariance(DIMENSION,'VariableNames',VNAMES) also 
%   specifies the variable names associated with rows and columns of this 
%   covariance matrix.
%
%   OBJ = IsotropicCovariance(DIMENSION,'VariableNames',VNAMES,'Name',NAME) 
%   also specifies the name to be associated with this covariance matrix.
        
        this = this@classreg.regr.lmeutils.covmats.CovarianceMatrix(dimension,varargin{:});
        
        this.DiagonalElements = logical(eye(this.Size));
        
        end % end of IsotropicCovariance
        
         function verifyTransformations(this)
%verifyTransformations - Check the parameter transformations between 
%   unconstrained, natural and canonical for accuracy.
          
            % (1) theta to D and D to theta
            theta = randn(1,1);
            D = theta2D(this,theta);            
            theta_recon = D2theta(this,D);            
            err(1) = max(abs(theta(:) - theta_recon(:)));            

            % (2) D to theta and theta to D
            D = rand(this.Size);
            D = mean(diag(D'*D))*eye(this.Size);
            theta = D2theta(this,D);
            D_recon = theta2D(this,theta);
            err(2) = max(abs(D(:) - D_recon(:)));
                                    
            % (3) 
            % (A) L via chol vs. getLowerTriangularCholeskyFactor
            theta = randn(1,1);
            D = theta2D(this,theta);            
            L = diag(sqrt(diag(D)));            
            this = setUnconstrainedParameters(this,theta);
            L1 = getLowerTriangularCholeskyFactor(this);            
            err(3) = max(abs(L(:) - L1(:)));

            % (3)
            % (B) theta to natural and then getLowerTriangularCholeskyFactor
            natural = convertUnconstrainedToNatural(this,theta);
            this = setNaturalParameters(this,natural);
            L2 = getLowerTriangularCholeskyFactor(this);            
            err(4) = max(abs(L(:) - L2(:)));
            
            % (3) 
            % (C) canonical parameters from unconstrained vs. natural
            this = setUnconstrainedParameters(this,theta);
            canonical1 = getCanonicalParameters(this);
            this = setNaturalParameters(this,natural);
            canonical2 = getCanonicalParameters(this);
            err(5) = max(abs(canonical1(:) - canonical2(:)));
                                    
            if max(err(:)) <= sqrt(eps)
                disp(getString(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:String_TransformsOK')));
            else
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:String_TransformsNotOK'));
            end            
            
        end % end of verifyTransformations.
        
    end

    methods (Access=public, Hidden=true)
        % A subclass must know how to initialize theta and how to convert
        % to L. For verifying transforms, it is also useful for subclasses
        % to define theta2D and D2theta methods.
        function theta = initializeTheta(this) %#ok<MANU>
        % Set D equal to eye(this.Size). This means theta would be a vector
        % of all zeros.    
            
            q = 1;
            theta = zeros(q, 1); 
        
        end % end of initializeTheta.
        
        function L = theta2L(this,theta) 
        
            L = sqrt(exp(theta))*eye(this.Size);
            
        end % end of theta2L.
        
        function D = theta2D(this,theta)
        
            D = exp(theta)*eye(this.Size);

        end % end of theta2D.
        
        function theta = D2theta(this,D) %#ok<INUSL>
        % We will assume that D is symmetric and positive definite. 
        
            theta = log(D(1,1));                    
            
        end % end of D2theta.
               
        % A subclass must know how to convert between unconstrained and
        % natural parameterizations.
        function natural = convertUnconstrainedToNatural(this,unconstrained)
%    Omega(i,i) = log( StdCorr(i,i) ) 
%               = log( sqrt(PSI(i,i)) )
%               = log( sqrt(sigma^2 * D(i,i)) )
%               = log( sigma * sqrt(D(i,i)) )
%               = log(sigma) + 1/2 * log(D(i,i))
%
%    natural = log(sigma) + (1/2)*unconstrained

           sigma = getSigma(this);
           
           natural = log(sigma) + 0.5*unconstrained;
            
        end % end of convertUnconstrainedToNatural.
        
        function unconstrained = convertNaturalToUnconstrained(this,natural)
%   PSI(i,i) = sigma^2 * D(i,i) and so
%    Omega(i,i) = log( StdCorr(i,i) ) 
%               = log( sqrt(PSI(i,i)) )
%               = log( sqrt(sigma^2 * D(i,i)) )
%               = log( sigma * sqrt(D(i,i)) )
%               = log(sigma) + 1/2 * log(D(i,i))
%     
%    log(D(i,i)) = 2*(Omega(i,i) - log(sigma))
%
%    unconstrained = 2*(natural - log(sigma))

            sigma = getSigma(this);
            
            unconstrained = 2*(natural - log(sigma));
 
        end % end of convertNaturalToUnconstrained.
        
        % A subclass must know how to convert between natural and canonical
        % parameterizations.
        function canonical = convertNaturalToCanonical(this,natural) %#ok<INUSL>
%   Suppose StdCorr is a matrix such that:
%              
%                StdCorr(i,i) = sqrt(PSI(i,i))
%                StdCorr(i,j) = 0
%
%   and suppose Omega is a matrix such that:
%
%               Omega(i,i) = log( StdCorr(i,i) )
%               Omega(i,j) = 0
%
%   then the Natural parameterization eta is given by diagonal elements of 
%   Omega and the Canonical parameterization heta is given by diagonal 
%   elements of StdCorr.

            canonical = exp(natural);      
            
        end % end of convertNaturalToCanonical.
        
        function natural = convertCanonicalToNatural(this,canonical) %#ok<INUSL>
%   Suppose StdCorr is a matrix such that:
%              
%                StdCorr(i,i) = sqrt(PSI(i,i))
%                StdCorr(i,j) = 0
%
%   and suppose Omega is a matrix such that:
%
%               Omega(i,i) = log( StdCorr(i,i) )
%               Omega(i,j) = 0
%
%   then the Natural parameterization eta is given by diagonal elements of 
%   Omega and the Canonical parameterization heta is given by diagonal 
%   elements of StdCorr.

            natural = log(canonical);            
            
        end % end of convertCanonicalToNatural.
        
        % A subclass must define canonical parameter names.
        function ds = defineCanonicalParameterNames(this)
        % Name the diagonal of StdCorr matrix. 
        %
        % If diagonal PSI is 3 by 3 with variable names {'v1','v2','v3'} 
        % then the elements of PSI can be named using the row and column 
        % variable names:
        %                  v1       v2       v3
        %     PSI = v1   psi11      0        0
        %           v2    0        psi22     0
        %           v3    0         0      psi33
        %
        % In this case, ds is a table array that looks like this:
        %
        %   Group      Name1  Name2  Type   Element of PSI
        %    Name        v1    v1    std         psi11        
        %    Name        v2    v2    std         psi22
        %    Name        v3    v3    std         psi33
        % 
        % where Name is the name of this Covariance matrix (this.Name).
        
            q = this.Size;                                    
            nrows = q*(q+1)/2;            
            group = cell(nrows,1);
            group(:) = {this.Name};            
            name1 = cell(nrows,1);
            name2 = cell(nrows,1);
            type = cell(nrows,1);  
                        
            vnames = this.VariableNames;
            count = 1;
            for col = 1:q
                for row = col:q                    
                    if (row == col)
                        name1{count} = vnames{row};
                        name2{count} = vnames{col};
                        type{count} = 'std';  
                        count = count + 1;
                    end
                end
            end
    
            ds = table(group,name1,name2,type,'VariableNames',{'Group','Name1','Name2','Type'});
            ds.Group = char(ds.Group);

            % Retain only the first row of ds since this is an isotropic
            % diagonal covariance.
            if q >= 1
                ds = ds(1,:);
            end
            
        end % end of defineCanonicalParameterNames.
        
        % A subclass must define its type (e.g., 'Full','Diagonal' etc.)
        function type = defineType(this) %#ok<MANU>
        % A 'Isotropic' covariance has type 'Isotropic'.
            type = classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_ISOTROPIC;
        end % end of defineType.
        
        % A subclass must define its covariance pattern.
        function covariancepattern = defineCovariancePattern(this)
        % A 'Isotropic' covariance sets the off diagonal elements to be 
        % equal to 0.         
            covariancepattern = logical(eye(this.Size));            
        end % end of defineCovariancePattern.        
        
    end
        
end

