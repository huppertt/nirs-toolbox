classdef CompoundSymmetryCovariance < classreg.regr.lmeutils.covmats.CovarianceMatrix
%  This feature is intended for internal use only and is subject to change
%  at any time without warning.

%CompoundSymmetryCovariance - Class to represent a 'Compound Symmetry' covariance matrix. 
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
%   When PSI is of Type 'Compound Symmetry', all diagonal elements of PSI 
%   are equal and all non-diagonal elements of PSI are equal. Because D is
%   a scaled version of PSI, D also has the 'Compound Symmetry' structure.
%   Here's an example to illustrate what the unconstrained, natural and 
%   canonical parameterizations look like in this case:
%
%   D has the following structure:
%
%   D = sigma1^2 * [1   rho rho ...
%                  rho   1  rho ...
%                  rho  rho  1  ...]
%
%   The smallest eigenvalue of D is sigma1^2*(1 +(q-1)*rho) where q is the
%   size of D. Positive semidefiniteness of D requires that:
%
%               -1/(q-1) <= rho <= 1   and
%                sigma1^2 >= 0
%   
%   ASIDE: If a <= x <= b then we can use the unconstrained parameter phi
%   like this:
%               x = a + exp(phi)/(1+exp(phi)) * (b-a)
%   As phi goes to -Inf, x goes to a and as phi goes to +Inf, x goes to b.
%   The inverse transform is:
%             phi = log( (x-a)/(b-x) )
%
%   To impose positive semidefinite compound symmetry on D, we use the
%   following parameterization:
%
%     sigma1^2 = exp(theta(1)) and
%          rho = -1/(q-1) + exp(theta(2))/(1+exp(theta(2)) * (1 + 1/(q-1))
%
%   In other words, the unconstrained parameters are:
%
%   theta(1) = log(sigma1^2)
%   theta(2) = log( (rho + 1/(q-1))/(1-rho) )
%
%   The canonical parameters are:
%
%   heta(1) = sqrt(sigma^2 * sigma1^2) and
%   heta(2) = rho
%
%   The natural parameters are:
%
%   eta(1) = log( heta(1) ) and
%   eta(2) = log( (1+heta(2))/(1-heta(2)) )
%
%   heta(1) = exp( eta(1) );
%   heta(2) = rho = (exp(eta(2)) - 1)/(exp(eta(2)) + 1)
%
%   CompoundSymmetryCovariance methods:
%       CompoundSymmetryCovariance - Create a 'Compound Symmetry' covariance matrix.
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
%   CompoundSymmetryCovariance properties:
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
        
         function this = CompoundSymmetryCovariance(dimension,varargin)
%CompoundSymmetryCovariance - Construct a CompoundSymmetryCovariance object.
%   OBJ = CompoundSymmetryCovariance(DIMENSION) takes an integer DIMENSION 
%   and constructs a IsotropicCovariance object.
%
%   OBJ = CompoundSymmetryCovariance(DIMENSION,'VariableNames',VNAMES) also 
%   specifies the variable names associated with rows and columns of this 
%   covariance matrix.
%
%   OBJ = CompoundSymmetryCovariance(DIMENSION,...,'Name',NAME) 
%   also specifies the name to be associated with this covariance matrix.
        
        this = this@classreg.regr.lmeutils.covmats.CovarianceMatrix(dimension,varargin{:});
        
        this.DiagonalElements = logical(eye(this.Size));
        
        end % end of IsotropicCovariance
        
         function verifyTransformations(this)
%verifyTransformations - Check the parameter transformations between 
%   unconstrained, natural and canonical for accuracy.
          
            % (1) theta to D and D to theta
            if this.Size == 1
                theta = randn(1,1);
            else
                theta = randn(2,1);
            end
            D = theta2D(this,theta);            
            theta_recon = D2theta(this,D);            
            err(1) = max(abs(theta(:) - theta_recon(:)));                        
            
            % (2) D to theta and theta to D
            q = this.Size;
            if q > 1
                rho = -1/(q-1) + 0.05;
                sigma1_sqrd = 3;
                D = sigma1_sqrd*(rho*ones(q) + (1-rho)*eye(q)); 
            else
                D = 3;
            end
            theta = D2theta(this,D);
            D_recon = theta2D(this,theta);
            err(2) = max(abs(D(:) - D_recon(:)));
                                    
            % (3) 
            % (A) L via chol vs. getLowerTriangularCholeskyFactor
            if this.Size == 1
                theta = randn;
            else
                theta = randn(2,1);
            end
            D = theta2D(this,theta);            
            L = chol(D,'lower');            
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
        function theta = initializeTheta(this)
        % Set D equal to eye(this.Size). This means theta would be a vector
        % of all zeros.    
            
            q = this.Size;
            if q == 1
                theta = zeros(1, 1); 
            else
                % sigma1_sqrd = 1 and rho = 0.
                theta = zeros(2,1);
                theta(1) = 0;
                theta(2) = -log(q-1);
            end
            
        end % end of initializeTheta.
        
        function L = theta2L(this,theta) 
        
            D = theta2D(this,theta); 
            try
                L = chol(D,'lower'); 
            catch ME %#ok<NASGU>
                L = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(D);                
            end 
            
        end % end of theta2L.
        
        function D = theta2D(this,theta)        
%     sigma1^2 = exp(theta(1)) and
%          rho = -1/(q-1) + exp(theta(2))/(1+exp(theta(2)) * (1 + 1/(q-1))
            
            if length(theta) == 1               
                sigma1_sqrd = exp(theta(1));
                D = sigma1_sqrd;                
            else
                sigma1_sqrd = exp(theta(1));
                q = this.Size;
                f = 1/(q-1);
                temp = exp(theta(2));
                rho = -f + (temp/(1+temp)) * (1 + f);
                
                D = sigma1_sqrd*(rho*ones(q) + (1-rho)*eye(q));                
            end
            
        end % end of theta2D.
        
        function theta = D2theta(this,D)
        % We will assume that D is symmetric and positive definite. 
        
            q = this.Size;            
            if q == 1
                sigma1_sqrd = D(1,1);
                theta(1) = log(sigma1_sqrd);
            else
                sigma1_sqrd = D(1,1);
                rho = D(1,2)/sigma1_sqrd;
            
                %   theta(1) = log(sigma1^2)
                %   theta(2) = log( (rho + 1/(q-1))/(1-rho) )
                       
                f = 1/(q-1);
            
                theta = zeros(2,1);
                theta(1) = log(sigma1_sqrd);
                theta(2) = log( (rho + f)/(1 - rho) );
            end            
            
        end % end of D2theta.
               
        % A subclass must know how to convert between unconstrained and
        % natural parameterizations.
        function natural = convertUnconstrainedToNatural(this,unconstrained)
%   sigma1^2 = exp(theta(1)) and
%   rho = -1/(q-1) + exp(theta(2))/(1+exp(theta(2)) * (1 + 1/(q-1))
%   heta(1) = sqrt(sigma^2 * sigma1^2) and
%   heta(2) = rho
%   eta(1) = log( heta(1) ) and
%   eta(2) = log( (1+heta(2))/(1-heta(2)) )

           sigma = getSigma(this);
           
           if length(unconstrained) == 1
               sigma1_sqrd = exp(unconstrained(1));
               canonical(1) = sqrt( (sigma^2) * sigma1_sqrd );
               natural(1) = log( canonical(1) );
           else
               q = this.Size;               
               f = 1/(q-1);
               temp = exp(unconstrained(2));
               rho = -f + (temp/(1+temp)) * (1 + f);
               natural = zeros(2,1);
               sigma1_sqrd = exp(unconstrained(1));
               canonical(1) = sqrt( (sigma^2) * sigma1_sqrd );
               natural(1) = log( canonical(1) );
               canonical(2) = rho;
               natural(2) = log( (1+canonical(2))/(1-canonical(2)) );
           end
            
        end % end of convertUnconstrainedToNatural.
        
        function unconstrained = convertNaturalToUnconstrained(this,natural)
%   heta(1) = exp( eta(1) );
%   heta(2) = rho = (exp(eta(2)) - 1)/(exp(eta(2)) + 1)
%   heta(1) = sqrt(sigma^2 * sigma1^2) and
%   heta(2) = rho
%   theta(1) = log(sigma1^2)
%   theta(2) = log( (rho + 1/(q-1))/(1-rho) )

            sigma = getSigma(this);
            canonical = convertNaturalToCanonical(this,natural);
            
            if length(canonical) == 1
               sigma1_sqrd = (canonical(1)^2)/(sigma^2);
               unconstrained(1) = log(sigma1_sqrd);               
            else
               sigma1_sqrd = (canonical(1)^2)/(sigma^2);
               rho = canonical(2);
               q = this.Size;
               unconstrained = zeros(2,1);
               unconstrained(1) = log(sigma1_sqrd); 
               unconstrained(2) = log( (rho + 1/(q-1))/(1-rho) );
            end
 
        end % end of convertNaturalToUnconstrained.
        
        % A subclass must know how to convert between natural and canonical
        % parameterizations.
        function canonical = convertNaturalToCanonical(this,natural) %#ok<INUSL>
%   heta(1) = exp( eta(1) );
%   heta(2) = rho = (exp(eta(2)) - 1)/(exp(eta(2)) + 1)

            if length(natural) == 1
                canonical(1) = exp(natural(1));
            else
                canonical = zeros(2,1);
                canonical(1) = exp(natural(1));
                
                % Split the next computation into 2 parts to avoid
                % overflow.
                if natural(2) <= 0
                    temp = exp(natural(2));
                    canonical(2) = (temp - 1)/(temp + 1);
                else
                    temp = exp(-natural(2));
                    canonical(2) = (1 - temp)/(1 + temp);
                end
            end
            
        end % end of convertNaturalToCanonical.
        
        function natural = convertCanonicalToNatural(this,canonical) %#ok<INUSL>
%   eta(1) = log( heta(1) ) and
%   eta(2) = log( (1+heta(2))/(1-heta(2)) )

            if length(canonical) == 1
                natural(1) = log(canonical(1));
            else
                natural = zeros(2,1);
                natural(1) = log(canonical(1));
                natural(2) = log( (1+canonical(2))/(1-canonical(2)) );
            end      
            
        end % end of convertCanonicalToNatural.
        
        % A subclass must define canonical parameter names.
        function ds = defineCanonicalParameterNames(this)
        % Name the elements of a 'Compound Symmetry' matrix.
        %
        % If PSI is 3 by 3 with variable names {'v1','v2','v3'} then the
        % elements of PSI can be named using the row and column variable
        % names:
        %                  v1       v2       v3
        %     PSI = v1   psi11    psi12     psi13
        %           v2   psi21    psi22     psi23
        %           v3   psi31    psi32     psi33
        %
        % In this case, ds is a table array that looks like this:
        %
        %   Group      Name1  Name2  Type   Element of PSI
        %    Name        v1    v1    std         psi11
        %    Name        v2    v1    corr        psi21
        %    Name        v3    v1    corr        psi31
        %    Name        v2    v2    std         psi22
        %    Name        v3    v2    corr        psi32
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
                    name1{count} = vnames{row};
                    name2{count} = vnames{col};
                    if (row == col)                                                
                        type{count} = 'std';                        
                    else                        
                        type{count} = 'corr';  
                    end
                    count = count + 1;
                end
            end
            
            ds = table(group,name1,name2,type,'VariableNames',{'Group','Name1','Name2','Type'});
            ds.Group = char(ds.Group);
            
            % If q = 1, retain only 1 row in ds otherwise retain only 2
            % rows since all other information is redundant and identical.
            if q == 1
                ds = ds(1,:);
            elseif q >= 2
                ds = ds(1:2,:);
            end
            
        end % end of defineCanonicalParameterNames.
        
        % A subclass must define its type (e.g., 'Full','Diagonal' etc.)
        function type = defineType(this) %#ok<MANU>
        % A 'Compound Symmetry' covariance has type 'Compound Symmetry'.
            type = classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_COMPSYMM;
        end % end of defineType.
        
        % A subclass must define its covariance pattern.
        function covariancepattern = defineCovariancePattern(this)
        % All elements of a 'Compound Symmetry' covariance matrix are
        % non-zero.
            covariancepattern = true(this.Size);            
        end % end of defineCovariancePattern.        
        
    end
        
end

