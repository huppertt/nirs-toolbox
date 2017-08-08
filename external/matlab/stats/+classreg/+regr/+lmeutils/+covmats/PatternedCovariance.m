classdef PatternedCovariance < classreg.regr.lmeutils.covmats.CovarianceMatrix
%  This feature is intended for internal use only and is subject to change
%  at any time without warning. 

%PatternedCovariance - Class to represent a 'Patterned' covariance matrix using the Cholesky parameterization.
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
%   When PSI is of Type 'Patterned', all elements of PSI are different
%   (except for symmetry) and some elements are constrained to be equal to
%   0. Here's an example to illustrate what the unconstrained, natural and
%   canonical parameterizations look like in this case:
%
%   Suppose D = [D_11 0    D_13
%                0    D_22 0
%                D_31 0    D_33]
%
%   The Cholesky factor of D is given by:
%
%           L = [L_11  0     0
%                0     L_22  0
%                L_31  0     L_33]
%
%   The Cholesky factor is not unique. An unconstrained parameterization
%   that guarantees positive definite D is given by the lower triangular
%   part of L excluding zeros in columnwise order.
%
%   In other words,
%
%           theta = [L_11;L_31;L_22;L_33] 
%
%   Suppose StdCorr is a matrix such that:
%              
%                StdCorr(i,i) = sqrt(PSI(i,i))
%                StdCorr(i,j) = PSI(i,j)/sqrt(PSI(i,i))/sqrt(PSI(j,j))
%
%   and suppose Omega is a matrix such that:
%
%               Omega(i,i) = log( StdCorr(i,i) )
%               Omega(i,j) = log( (1 + StdCorr(i,j))/(1 - StdCorr(i,j)) )
%
%   then the Natural parameterization eta is given by elements in the lower
%   triangular part of Omega excluding zeros and the Canonical
%   parameterization heta is given by elements in the lower triangular part
%   of StdCorr excluding zeros.
%
%   PatternedCovariance methods:
%       PatternedCovariance - Create a 'Patterned' covariance matrix.
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
%   PatternedCovariance properties:
%       Name - Name of this covariance matrix
%       Size - Size of this covariance matrix
%       VariableNames - Names of random variables whose covariance is PSI
%       Type - The type of this covariance matrix (e.g., "Diagonal").
%       CovariancePattern - Pattern of non-zeros in this covariance matrix
%       WhichParameterization - Parameterization that is currently in use.
%       NumParametersExcludingSigma - Number of parameters excluding sigma.    
%
%   See also classreg.regr.lmeutils.covmats.CovarianceMatrix
    
%   Copyright 2012-2014 The MathWorks, Inc.

    properties (Access=protected)
        
%An indicator matrix for the diagonal elements. Basically this is equal to
%   logical(eye(this.Size)).        
        DiagonalElements
        
    end
    
    methods (Access=public)
       
        function this = PatternedCovariance(dimension,varargin)
%PatternedCovariance - Construct a PatternedCovariance object.
%   OBJ = PatternedCovariance(DIMENSION) takes an integer DIMENSION and 
%   constructs a PatternedCovariance object.
%
%   OBJ = PatternedCovariance(DIMENSION,'VariableNames',VNAMES) also
%   specifies the variable names associated with rows and columns of this
%   covariance matrix.
%
%   OBJ = PatternedCovariance(DIMENSION,'VariableNames',VNAMES,'Name',NAME) 
%   also specifies the name to be associated with this covariance matrix.
%
%   OBJ = PatternedCovariance(DIMENSION,...,'CovariancePattern',PAT) also 
%   specifies a logical square matrix of size DIMENSION that specifies the
%   pattern of non-zeros in the PatternCovariance matrix. The default value
%   for 'CovariancePattern' is an all true matrix of size DIMENSION.
        
        this = this@classreg.regr.lmeutils.covmats.CovarianceMatrix(dimension,varargin{:});        
        this.DiagonalElements = logical(eye(this.Size));
                        
         % Parse CovariancePattern.
            dfltCovariancePattern = [];
            parnames = {'CovariancePattern'};
            dflts = {dfltCovariancePattern};
            [covariancepattern,~,~] = internal.stats.parseArgs(parnames,dflts,varargin{:});

            if isempty(covariancepattern)
                covariancepattern = true(this.Size);
            else
                msg = ['''CovariancePattern'' must be a logical matrix of size ',num2str(this.Size)];
                
                assert( ismatrix(covariancepattern) & ...
                    islogical(covariancepattern) & ...
                    all(size(covariancepattern) == [this.Size,this.Size]), msg );               
            end

            % Set Type and CovariancePattern.
            this.Type = defineType(this);      
            this.CovariancePattern = covariancepattern;
            
            % Set initial sigma and theta.            
            theta = initializeTheta(this);
            this.Sigma = 1; 
            
            % Set unconstrained parameters and parameterization type
            this.UnconstrainedParameters = theta;
            this.WhichParameterization = classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_UNCONSTRAINED;
            
            % How many parameters exclusing sigma? What are canonical
            % parameter names?
            this.NumParametersExcludingSigma = length(theta);
            this.CanonicalParameterNames = defineCanonicalParameterNames(this);
            
            % Ensure that each unknown parameter has a name and that the
            % number of names matches the number of unknown parameters.
            assert( size(this.CanonicalParameterNames,1) == length(theta) );                   
            
        end % end of PatternedCovariance
        
        function verifyTransformations(this)
%verifyTransformations - Check the parameter transformations between 
%   unconstrained, natural and canonical for accuracy.
          
            % (1) theta to D and D to theta
            theta = generateRandomTheta(this);
            D = theta2D(this,theta);            
            theta_recon = D2theta(this,D);            
            err(1) = max(abs(theta(:) - theta_recon(:)));            

            % (2) D to theta and theta to D
            D = rand(this.Size);
            D = D'*D;
            D(~this.CovariancePattern) = 0;
            theta = D2theta(this,D);
            D_recon = theta2D(this,theta);
            err(2) = max(abs(D(:) - D_recon(:)));
                                    
            % (3) 
            % (A) L via chol vs. getLowerTriangularCholeskyFactor
            theta = generateRandomTheta(this);
            D = theta2D(this,theta);            
            L = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(D);            
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
        % Set D equal to eye(this.Size). This means L would also be 
        % eye(this.Size) and so theta would be the lower triangular part of
        % the identity matrix eye(this.Size).    
            
            q = this.Size;
            L = eye(q);
            theta = L(tril(this.CovariancePattern));            
        
        end % end of initializeTheta.
        
        function L = theta2L(this,theta)
        
            L = zeros(this.Size);
            L(tril(this.CovariancePattern)) = theta;           
            
        end % end of theta2L.
        
        function D = theta2D(this,theta)
        
            L = theta2L(this,theta);
            D = L*L';
            
        end % end of theta2D.
        
        function theta = D2theta(this,D)
        % We will assume that D is symmetric and positive definite. 
        
            try
               % Use standard chol.
               L = chol(D,'lower');                               
            catch ME %#ok<NASGU>
               % Use singularLowerChol.
               L = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(D);
            end            
            % Extract lower triangular part of L into theta.
            theta = L(tril(this.CovariancePattern));
            
        end % end of D2theta.
               
        % A subclass must know how to convert between unconstrained and
        % natural parameterizations.
        function natural = convertUnconstrainedToNatural(this,unconstrained)
            
           % First get PSI.
           sigma = getSigma(this);
           D = theta2D(this,unconstrained);
           PSI = (sigma^2)*D;
           
           % Next get StdCorr followed by Omega.
           StdCorr = PSI2StdCorr(this,PSI);
           Omega = StdCorr2Omega(this,StdCorr);
           
           % Extract the lower triangular part of Omega.
           natural = Omega(tril(this.CovariancePattern));           
            
        end % end of convertUnconstrainedToNatural.
        
        function unconstrained = convertNaturalToUnconstrained(this,natural)
           
            % First make Omega.
            Omega = zeros(this.Size);
            Omega(tril(this.CovariancePattern)) = natural;
            % Make Omega symmetric.
            Omega = Omega + tril(Omega,-1)';
            
            % Convert Omega to StdCorr.
            StdCorr = Omega2StdCorr(this,Omega);
            % Convert StdCorr to PSI.
            PSI = StdCorr2PSI(this,StdCorr);
            
            % Get D from PSI.
            sigma = getSigma(this);
            D = PSI/(sigma^2);
            
            % Finally, get theta from D.
            unconstrained = D2theta(this,D);
 
        end % end of convertNaturalToUnconstrained.
        
        % A subclass must know how to convert between natural and canonical
        % parameterizations.
        function canonical = convertNaturalToCanonical(this,natural)
           
            % First make Omega.
            Omega = zeros(this.Size);
            Omega(tril(this.CovariancePattern)) = natural;
            % Make Omega symmetric.
            Omega = Omega + tril(Omega,-1)';
            
            % Convert Omega to StdCorr.
            StdCorr = Omega2StdCorr(this,Omega);
            
            % Extract the lower triangular part of StdCorr.            
            canonical = StdCorr(tril(this.CovariancePattern));            
            
        end % end of convertNaturalToCanonical.
        
        function natural = convertCanonicalToNatural(this,canonical)
            
            % First make StdCorr.
            StdCorr = zeros(this.Size);
            StdCorr(tril(this.CovariancePattern)) = canonical;
            % Make StdCorr symmetric.
            StdCorr = StdCorr + tril(StdCorr,-1)';
            
            % Convert StdCorr to Omega.
            Omega = StdCorr2Omega(this,StdCorr);
            
            % Extract the lower triangular part of Omega.
            natural = Omega(tril(this.CovariancePattern));
            
        end % end of convertCanonicalToNatural.
        
        % A subclass must define canonical parameter names.
        function ds = defineCanonicalParameterNames(this)
        % Name the lower triangular part of StdCorr matrix. 
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
        
            group = cell(this.NumParametersExcludingSigma,1);
            group(:) = {this.Name};
            
            name1 = cell(this.NumParametersExcludingSigma,1);
            name2 = cell(this.NumParametersExcludingSigma,1);
            type = cell(this.NumParametersExcludingSigma,1);  
                        
            q = this.Size;
            vnames = this.VariableNames;
            count = 1;
            for col = 1:q
                for row = col:q
                    if this.CovariancePattern(row,col) == true
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
            end            
            
            ds = table(group,name1,name2,type,'VariableNames',{'Group','Name1','Name2','Type'});
            ds.Group = char(ds.Group);

        end % end of defineCanonicalParameterNames.
        
        % A subclass must define its type (e.g., 'FullCholesky','Diagonal' etc.)
        function type = defineType(this) %#ok<MANU>
        % A 'Patterned' covariance has type 'Patterned'.
            type = classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_PATTERNED;
        end % end of defineType.
        
        % A subclass must define its covariance pattern.
        function covariancepattern = defineCovariancePattern(this)
            if isempty(this.CovariancePattern)
                covariancepattern = true(this.Size);
            else
                covariancepattern = this.CovariancePattern;  
            end
        end % end of defineCovariancePattern.
        
    end
    
    methods (Access=private)           
        
        function theta = generateRandomTheta(this,condnum)
        % A utility function to generate a sensible random value of the 
        % unconstrained parameter vector. You can specify the desired 
        % condition number of D in condnum. The default value is 1e5.
            
            if nargin < 2
                condnum = 1e5;
            end        
            % Generate elements of the lower triangular Cholesky factor
            A = randn(this.Size,this.Size);
            Q = orth(A);
            lambda = linspace(1,condnum,this.Size);
            D = Q*diag(lambda)*Q';
            theta = D2theta(this,D);                        

        end % end of generateRandomTheta.
        
        function StdCorr = PSI2StdCorr(this,PSI)
% StdCorr(i,i) = sqrt(PSI(i,i))
% StdCorr(i,j) = PSI(i,j)/sqrt(PSI(i,i))/sqrt(PSI(j,j))
  
            [StdCorr,SIGMA] = corrcov(PSI);
            StdCorr(this.DiagonalElements) = SIGMA;                       

        end % end of PSI2StdCorr.

        function Omega = StdCorr2Omega(this,StdCorr)
% Omega(i,i) = log( StdCorr(i,i) )
% Omega(i,j) = log( (1 + StdCorr(i,j))/(1 - StdCorr(i,j)) )
            
            Omega = log( (1+StdCorr)./(1-StdCorr) );
            diagelem = this.DiagonalElements;
            Omega(diagelem) = log( StdCorr(diagelem) );
            
        end % end of StdCorr2Omega.

        function StdCorr = Omega2StdCorr(this,Omega)
% StdCorr(i,i) = exp( Omega(i,i) )
% StdCorr(i,j) = ( exp(Omega(i,j)) - 1 )/( exp(Omega(i,j)) + 1 )

            % For non-diagonal elements.
            temp = exp(Omega);
            StdCorr = (temp - 1)./(temp + 1);
            
            % For diagonal elements.
            diagelem = this.DiagonalElements;
            StdCorr(diagelem) = exp(Omega(diagelem)); 

        end % end of Omega2StdCorr.

        function PSI = StdCorr2PSI(this,StdCorr)
% StdCorr(i,i) = sqrt(PSI(i,i))
% StdCorr(i,j) = PSI(i,j)/sqrt(PSI(i,i))/sqrt(PSI(j,j))

            SIGMA = diag(StdCorr);
            StdCorr(this.DiagonalElements) = 1;
            
            % Multiply column j by SIGMA(j)
            StdCorr = bsxfun(@times,StdCorr,SIGMA');
            
            % Multiply row i by SIGMA(i)
            PSI = bsxfun(@times,StdCorr,SIGMA);

        end % end of StdCorr2PSI.                
        
    end
    
end

