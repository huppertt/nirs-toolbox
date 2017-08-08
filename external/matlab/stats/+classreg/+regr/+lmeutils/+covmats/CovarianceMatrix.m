classdef CovarianceMatrix
%  This feature is intended for internal use only and is subject to change
%  at any time without warning. 

%CovarianceMatrix - Abstract class for parameterized covariance matrices.
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
%   (A) During LinearMixedModel fitting:
%        
%       (1) We normally profile out sigma^2 and so the ML or REML objective 
%           function depends only on theta through L(theta).
%       (2) To initialize optimization, we need an initial value for theta,
%           say theta0.
%       (3) To perform optimization, we need to transform theta into
%           L(theta) repeatedly.
%
%   (B) After optimization completes, we get the ML or REML estimate of
%   theta, say theta_hat. Using this, we also get the ML or REML estimate
%   of sigma^2, say sigma^2_hat. We would like to compute the covariance of
%   [theta_hat; log(sigma_hat)]. For this, we express the log likelihood 
%   profiled on fixed effects (for ML) or the restricted log likelihood 
%   (for REML) as a function of [theta;log(sigma)]. Inside this function, 
%   we need to convert from [theta;log(sigma)] to L(theta) and sigma^2.
%
%   (C) To compute confidence intervals (CIs) on variance components, we
%   re-express the log likelihood profiled on the fixed effects (for ML) or
%   the restricted log likelihood as a function of [eta;log(sigma)] where
%   eta is a "Natural" parameter vector. The parameters in eta are
%   individually unconstrained but not jointly so in such a way that eta
%   can be transformed into sensible "Canonical" parameters by 
%   monotonically increasing functions h_i(eta_i). When deriving CIs for
%   [eta;log(sigma)], we have to convert from [eta;log(sigma)] to L(theta)
%   and sigma^2.
%
%   Suppose we derive CIs on eta_i and log(sigma) that look like this:
%
%               a <= eta_i <= b
%               c <= log(sigma) <= d
%
%   Applying the monotonically increasing transformation h_i(eta_i) we get:
%
%             h_i(a) <= h_i(eta_i) <= h_i(b)
%             exp(c) <= sigma <= exp(d)
%   
%   The h_i(eta_i) are the "Canonical" parameters. In other words, we need
%   the transformation [eta;log(sigma)] to [h(eta);sigma] where h(eta) =
%   [h_1(eta_1);h_2(eta_2);...].
%
%   Each covariance matrix has 3 parameterizations:
% 
%   (1) The unconstrained parameterization: This is formed from the
%       unconstrained parameters theta and log(sigma).
%
%   (2) The Natural parameterization: This is formed from the Natural
%       parameters eta and log(sigma).
%
%   (3) The Canonical parameterization. This is formed from the canonical
%       parameters h(eta) and sigma where elements h_i are monotonically 
%       increasing functions of eta_i and h_i(eta_i) "make sense" for 
%       interpretation purposes.
%
%   This class provides an interface for writing various parameterized
%   covariance matrices and provides functions for:
%
%   (1) constructing a covariance matrix.
%   (2) getting unconstrained, natural and canonical parameterizations.
%   (3) setting unconstrained and natural parameterizations.
%   (4) converting between Natural and unconstrained parameterization.
%   (5) converting between Natural and Canonical parameterization.
%   (6) converting between Canonical and unconstrained parameterization.
%   (7) getting/setting sigma.
%   (8) getting L(theta).
%
%   CovarianceMatrix methods:
%       CovarianceMatrix - Create a CovarianceMatrix object.
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
%   CovarianceMatrix properties:
%       Name - Name of this covariance matrix
%       Size - Size of this covariance matrix
%       VariableNames - Names of random variables whose covariance is PSI
%       Type - The type of this covariance matrix (e.g., "Diagonal")
%       CovariancePattern - Pattern of non-zeros in this covariance matrix
%       WhichParameterization - Parameterization that is currently in use
%       NumParametersExcludingSigma - Number of parameters excluding sigma.    

%   Copyright 2012-2014 The MathWorks, Inc.

    properties (Access=public)
%Name - A string (a tag) associated with this covariance matrix.
        Name        
    end

    properties (GetAccess=public, SetAccess=protected)        
%Size - The number of rows and columns in this covariance matrix.        
        Size
        
%VariableNames - A cell array containing the names of the random variables
%   whose covariance matrix is represented by this class.
        VariableNames        
        
%Type - The type of this covariance matrix. For example if Type = 'Full'
%   then we have a full covariance matrix. If Type = 'Diagonal' then we
%   have a diagonal covariance matrix.
        Type
        
%CovariancePattern - A logical matrix indicating the pattern of non-zeros
%   in this covariance matrix. If CovariancePattern(a,b) = false then
%   element (a,b) of the covariance matrix is set to 0.
        CovariancePattern       
        
%WhichParameterization - The parameterization that is currently in use.
        WhichParameterization  
        
%NumParametersExcludingSigma - Number of parameters excluding sigma.    
        NumParametersExcludingSigma                      
    end
    
    properties (Access=protected)
%UnconstrainedParameters - The unconstrained parameter vector associated 
%   with this covariance matrix.      
        UnconstrainedParameters
                  
%NaturalParameters - The Natural parameter vector associated with this
%   covariance matrix.              
        NaturalParameters
    
%CanonicalParameters - The Canonical parameter vector associated with this
%   covariance matrix.
        CanonicalParameters        
        
%CanonicalParameterNames - A cell array containing the name associated with
%   each element of CanonicalParameters.    
        CanonicalParameterNames                   
        
%Sigma - The sigma value associated with this covariance matrix.
        Sigma
 
%SigmaName - Name of the sigma parameter.
        SigmaName = 'Res std';
    end
    
    properties (Access=public,Constant=true,Hidden=true)      
        % Various parameterizations that we support.
        PARAMETERIZATION_NATURAL = 'Natural';
        PARAMETERIZATION_CANONICAL = 'Canonical';
        PARAMETERIZATION_UNCONSTRAINED = 'Unconstrained';    
        
        % Various covariance matrix structures that we support.
        TYPE_FULL = 'Full';
        TYPE_FULLCHOLESKY = 'FullCholesky';
        TYPE_DIAGONAL = 'Diagonal';
        TYPE_ISOTROPIC = 'Isotropic';
        TYPE_COMPSYMM = 'CompSymm';
        TYPE_BLOCKED = 'Blocked';
        TYPE_FIXEDWEIGHTS = 'FixedWeights';
        TYPE_PATTERNED = 'Patterned';
        
        % Allowed covariance structures. To create a BlockedCovariance, use
        % the BlockedCovariance constructor. Otherwise, you can use the
        % createCovariance convenience method.
        AllowedCovarianceTypes = ...
           {classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FULL,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FULLCHOLESKY,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_DIAGONAL,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_ISOTROPIC,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_COMPSYMM,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FIXEDWEIGHTS,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_PATTERNED};
    end
    
    methods (Access=public)        
                
        function this = CovarianceMatrix(dimension,varargin)
%CovarianceMatrix - Construct a CovarianceMatrix object.
%   OBJ = CovarianceMatrix(DIMENSION) takes an integer DIMENSION and 
%   constructs a CovarianceMatrix object.
%
%   OBJ = CovarianceMatrix(DIMENSION,'VariableNames',VNAMES) also specifies 
%   the variable names associated with rows and columns of this covariance
%   matrix.
%
%   OBJ = CovarianceMatrix(DIMENSION,'VariableNames',VNAMES,'Name',NAME) 
%   also specifies the name to be associated with this covariance matrix.


            if nargin == 0
                return;
            end

            % Provide default values for optional parameters.
            dfltVariableNames = [];
            dfltName = [];

            % Set defaults for optional parameters.
            parnames = {'VariableNames'  ,   'Name'  };
               dflts = {dfltVariableNames,   dfltName};
                            
            % Process optional param name/value pairs.
            [variablenames,name,~,~] = internal.stats.parseArgs(parnames,dflts,varargin{:});
            
            % Validate dimension.
            isPosInteger = internal.stats.isIntegerVals(dimension,0);
            if ~isPosInteger || ~isscalar(dimension)
                % <entry key="BadDimension">DIMENSION must be a scalar positive integer.</entry>
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadDimension'));
            end
            
            % If supplied, validate variablenames.
            if ~isempty(variablenames)
                % Make sure variablenames is a cell vector of strings.
                if ~internal.stats.isStrings(variablenames,true) || ~isvector(variablenames)
                    % <entry key="BadVariableNames">''VariableNames'' must be specified as a cell vector of strings.</entry>
                    error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadVariableNames'));
                end

                % Make variablenames into a column vector
                if size(variablenames,1) == 1
                    variablenames = variablenames';
                end

                % Make sure columnnames have length equal to dimension
                if length(variablenames) ~= dimension
                    % <entry key="BadVariableNamesLength">''VariableNames'' must be of length {0}.</entry>
                    error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadVariableNamesLength',num2str(dimension)));
                end
            else
               % variablenames is empty. Assign default variable names.
               variablenames = cell(dimension,1);
               for i = 1:dimension
                   variablenames{i} = ['v',num2str(i)];
               end
            end
            
            % If supplied, validate name.
            if ~isempty(name)               
                % Make sure name is a string
                [tf,name] = internal.stats.isString(name,true);                
                if ~tf
                    % <entry key="BadName">''Name'' must be a string.</entry>
                    error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadName'));
                end                
            else
               % name is empty. Assign default name.
               name = 'g';
            end
            
            this.Name = name;
            this.Size = dimension;
            this.VariableNames = variablenames;                        
            
            this.Type = defineType(this);      
            assert( internal.stats.isString(this.Type,false) == true );
            
            this.CovariancePattern = defineCovariancePattern(this);
            assert( all(size(this.CovariancePattern) == [dimension,dimension]) );            
            
            % Set initial sigma and theta.
            this.WhichParameterization = classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_UNCONSTRAINED;
            theta = initializeTheta(this);
            this.UnconstrainedParameters = theta;
            
            this.NumParametersExcludingSigma = length(theta);

            this.CanonicalParameterNames = defineCanonicalParameterNames(this);
            
            % Ensure that each unknown parameter has a name and that the
            % number of names matches the number of unknown parameters.
            if ~isa(this,'classreg.regr.lmeutils.covmats.BlockedCovariance')
                assert( size(this.CanonicalParameterNames,1) == length(theta) );
            end
            
            this.Sigma = 1;                                    
            
        end % end of CovarianceMatrix.                        
        
        function unconstrained = getUnconstrainedParameters(this)
%getUnconstrainedParameters - Get unconstrained parameterization.
%   UNCONSTRAINED = getUnconstrainedParameters(OBJ) accepts an object OBJ
%   of type CovarianceMatrix and returns its unconstrained parameters.
            
            % Create unconstrained parameterization if required.
            if isempty(this.UnconstrainedParameters)
                this = createUnconstrainedParameterization(this);
            end            
            unconstrained = this.UnconstrainedParameters;           
            
        end % end of getUnconstrainedParameters.
        
        function natural = getNaturalParameters(this)
%getNaturalParameters - Get Natural parameterization.
%   NATURAL = getNaturalParameters(OBJ) accepts an object OBJ of type
%   CovarianceMatrix and returns its Natural parameters.
                        
            % Create natural parameterization if required.
            if isempty(this.NaturalParameters)
                this = createNaturalParameterization(this);
            end                       
            natural = this.NaturalParameters;
            
        end % end of getNaturalParameters.
        
        function canonical = getCanonicalParameters(this)
%getCanonicalParameters - Get Canonical parameterization.
%   CANONICAL = getCanonicalParameters(OBJ) accepts an object OBJ of type
%   CovarianceMatrix and returns its Canonical parameters.
                        
            % Create canonical parameterization if required.
            if isempty(this.CanonicalParameters)
                this = createCanonicalParameterization(this);
            end
            canonical = this.CanonicalParameters;
                                    
        end % end of getCanonicalParameters.
        
        function names = getCanonicalParameterNames(this)
%getCanonicalParameterNames - Get Canonical parameter names.
%   NAMES = getCanonicalParameterNames(OBJ) accepts an object OBJ of type
%   CovarianceMatrix and returns a cell array of strings or a table array
%   containing information about the Canonical parameter names.
                                   
            names = this.CanonicalParameterNames;       
            
        end % end of getCanonicalParameterNames.        
        
        function this = setUnconstrainedParameters(this,theta)
%setUnconstrainedParameters - Sets the unconstrained parameterization.
%   OBJ = setUnconstrainedParameters(OBJ,THETA) accepts an object OBJ of 
%   type CovarianceMatrix, an unconstrained parameter vector THETA and 
%   returns the object OBJ with the specified unconstrained parameters. 
%   The length of THETA must match the number of unconstrained parameters
%   (excluding sigma) in OBJ.

            narginchk(2,Inf);
            
            % Make sure theta is a real vector of the right size. 
            isRealVector = isreal(theta) & isvector(theta);            
            if isRealVector && length(theta) == this.NumParametersExcludingSigma             
                % Make theta into a column vector.
                if size(theta,1) == 1
                    theta = theta';
                end                                                 
            else
                % <entry key="BadTheta">theta must be a real vector of length {0}.</entry>
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadTheta',num2str(this.NumParametersExcludingSigma)));
            end

            % Forget Natural and Canonical parameterizations.
            this.NaturalParameters = [];
            this.CanonicalParameters = [];
            
            % Set the unconstrained parameterization.
            this.UnconstrainedParameters = theta;            
            this.WhichParameterization = classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_UNCONSTRAINED;
            
        end % end of setUnconstrainedParameters.
                       
        function this = setNaturalParameters(this,eta)
%setNaturalParameters - Sets the Natural parameterization.
%   OBJ = setNaturalParameters(OBJ,ETA) accepts an object OBJ of type
%   CovarianceMatrix, a Natural parameter vector ETA and returns the object
%   OBJ with the specified Natural parameters. The length of ETA must match
%   the number of unconstrained parameters (excluding sigma) in OBJ.

            narginchk(2,Inf);

            % Make sure eta is a real vector of the right size. 
            isRealVector = isreal(eta) & isvector(eta);            
            if isRealVector && length(eta) == this.NumParametersExcludingSigma             
                % Make theta into a column vector.
                if size(eta,1) == 1
                    eta = eta';
                end                                                 
            else
                % <entry key="BadEta">eta must be a real vector of length {0}.</entry>
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadEta',num2str(this.NumParametersExcludingSigma)));
            end
           
            % Forget Unconstrained and Canonical parameterizations.
            this.UnconstrainedParameters = [];
            this.CanonicalParameters = [];
            
            % Set the Natural parameterization.
            this.NaturalParameters = eta;           
            this.WhichParameterization = classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_NATURAL;
                        
        end % end of setNaturalParameters.        
        
        function L = getLowerTriangularCholeskyFactor(this)
%getLowerTriangularCholeskyFactor - Gets lower triangular Cholesky factor.
%   L = getLowerTriangularCholeskyFactor(OBJ) accepts an object of type
%   CovarianceMatrix and extracts its lower triangular Cholesky factor. If
%   PSI = sigma^2 * D is the covariance matrix represented by OBJ, then the
%   factor L is such that L*L' = D.

            narginchk(1,Inf);            

            % Get theta and transform to L.
            theta = getUnconstrainedParameters(this);                                                
                L = theta2L(this,theta);            
            
        end % end of getLowerTriangularCholeskyFactor.                                   
        
        function sigma = getSigma(this)    
%getSigma - Get sigma for this covariance matrix.
%   sigma = getSigma(OBJ) accepts an object of type CovarianceMatrix and 
%   extracts the sigma for this covariance matrix.

            sigma = this.Sigma;            

        end % end of getSigma.                
        
        function this = setSigma(this,sigma)
%setSigma - Set sigma for this covariance matrix.
%   OBJ = setSigma(OBJ,SIGMA) accepts an object of type CovarianceMatrix, 
%   a real positive scalar SIGMA, sets the Sigma property of OBJ and 
%   returns the modified object.

            narginchk(2,Inf);
            
            % Ensure sigma is a real positive scalar
            if ~isreal(sigma) || ~isscalar(sigma)
                % <entry key="BadSigma">SIGMA must be a real positive scalar.</entry>
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadSigma'));
            end
                        
            this.Sigma = sigma;            
            
        end % end of setSigma.           
       
        function this = createUnconstrainedParameterization(this)
%createUnconstrainedParameterization - Create unconstrained parameterization
%   OBJ = createUnconstrainedParameterization(OBJ) takes an object OBJ of
%   type CovarianceMatrix, stores the unconstrained parameterization in it
%   and returns the modified object.

            switch this.WhichParameterization   
                
                    case classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_NATURAL
                        
                        natural = this.NaturalParameters;
                        unconstrained = convertNaturalToUnconstrained(this,natural);
                        this.UnconstrainedParameters = unconstrained;
                        
                    case classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_CANONICAL
                        
                        canonical = this.CanonicalParameters;
                        unconstrained = convertCanonicalToUnconstrained(this,canonical);
                        this.UnconstrainedParameters = unconstrained;
                                                
                    case classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_UNCONSTRAINED
                        % Already using the parameterization
                    otherwise
                        % <entry key="BadParameterization">Only Unconstrained, Natural and Canonical parameterizations are supported.</entry>      
                        error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadParameterization'));
            end                            
        end % end of createUnconstrainedParameterization.
        
        function this = createNaturalParameterization(this)
%createNaturalParameterization - Create Natural parameterization.
%   OBJ = createNaturalParameterization(OBJ) takes an object OBJ of type
%   CovarianceMatrix, stores the Natural parameterization in it and returns
%   the modified object.
            
            switch this.WhichParameterization 
                
                    case classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_NATURAL
                        % Already using the parameterization
                    case classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_CANONICAL
                        
                        canonical = this.CanonicalParameters;
                        natural = convertCanonicalToNatural(this,canonical);
                        this.NaturalParameters = natural;
                        
                    case classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_UNCONSTRAINED
                        
                        unconstrained = this.UnconstrainedParameters;
                        natural = convertUnconstrainedToNatural(this,unconstrained);
                        this.NaturalParameters = natural;
                                               
                otherwise
                        % <entry key="BadParameterization">Only Unconstrained, Natural and Canonical parameterizations are supported.</entry>      
                        error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadParameterization'));                        
            end                    
        end % end of createNaturalParameterization.
        
        function this = createCanonicalParameterization(this)
%createCanonicalParameterization - Create Canonical parameterization.
%   OBJ = createCanonicalParameterization(OBJ) takes an object OBJ of type
%   CovarianceMatrix, stores the Canonical parameterization in it and
%   returns the modified object.

            switch this.WhichParameterization                  
                
                case classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_NATURAL
                    
                    natural = this.NaturalParameters;
                    canonical = convertNaturalToCanonical(this,natural);
                    this.CanonicalParameters = canonical;
                                            
                case classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_CANONICAL
                    % Already using the parameterization
                case classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_UNCONSTRAINED
                    
                    unconstrained = this.UnconstrainedParameters;
                    canonical = convertUnconstrainedToCanonical(this,unconstrained);
                    this.CanonicalParameters = canonical;
                                            
                otherwise
                    % <entry key="BadParameterization">Only Unconstrained, Natural and Canonical parameterizations are supported.</entry>      
                    error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadParameterization'));                    
            end            
        end % end of createCanonicalParameterization.               
        
    end
    
    methods (Access=protected)
        
        function canonical = convertUnconstrainedToCanonical(this,unconstrained)
%convertUnconstrainedToCanonical - convert unconstrained parameterization to canonical.           
            
            natural = convertUnconstrainedToNatural(this,unconstrained);
            canonical = convertNaturalToCanonical(this,natural);

        end % end of convertUnconstrainedToCanonical.
                
        function unconstrained = convertCanonicalToUnconstrained(this,canonical)
%convertCanonicalToUnconstrained - convert Canonical parameterization to unconstrained.

            natural = convertCanonicalToNatural(this,canonical);                        
            unconstrained = convertNaturalToUnconstrained(this,natural);        

        end % end of convertCanonicalToUnconstrained.
        
    end
    
    methods (Abstract=true, Access=public, Hidden=true)
      
        % A subclass must know how to initialize theta and how to convert
        % to L. For verifying transforms, it is also useful for subclasses
        % to define theta2D and D2theta methods.
        theta = initializeTheta(this);        
        L = theta2L(this,theta);        
        D = theta2D(this,theta);
        theta = D2theta(this,D);
               
        % A subclass must know how to convert between unconstrained and
        % natural parameterizations.
        natural = convertUnconstrainedToNatural(this,unconstrained);
        unconstrained = convertNaturalToUnconstrained(this,natural); 
        
        % A subclass must know how to convert between natural and canonical
        % parameterizations.
        canonical = convertNaturalToCanonical(this,natural);
        natural = convertCanonicalToNatural(this,canonical);                        
        
        % A subclass must define canonical parameter names.
        names = defineCanonicalParameterNames(this);       
        
        % A subclass must define its type (e.g., 'Full','Diagonal' etc.)
        type = defineType(this);            
        
        % A subclass must define its covariance pattern.
        covariancepattern = defineCovariancePattern(this);
        
    end
    
    methods (Access=public, Static=true)
        
        %=== helper function singularLowerChol
        function R1 = singularLowerChol(T)
%singularLowerChol - Cholesky for symmetric positive semidefinite matrix.
%   R1 = singularLowerChol(T) takes T, a real square matrix that is known
%   to be symmetric and positive semidefinite and computes its lower
%   triangular Cholesky factor R1. Here's how singularLowerChol works:
%
%   (1) T may not be exactly symmetric due to numerical roundoff, so T is 
%   first replaced by (T+T')/2.
% 
%   (2) If T contains NaNs or Infs then R1 is an all NaN lower triangular
%   matrix.
%
%   (3) The smallest eigenvalue of T may be a small negative number. It is
%   assumed that this is due to numerical roundoff. R1 is then taken to be
%   lower triangular Cholesky factor of (T + delta*I) where I is the
%   identity matrix and delta is a number such that the smallest eigenvalue
%   of (T + delta*I) is positive.
       
        % (1) Ensure that T is a real and square matrix.        
        % <entry key="BadT">T must be a real and square matrix.</entry>      
        assertThat(isreal(T) & ismatrix(T) & size(T,1) == size(T,2),'stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadT');
        
        % (2) Ensure that T is symmetric.
        T = (T+T')/2;
        
        % (3) Ensure that T has no NaNs and Infs. If yes, just return a
        % lower triangular matrix of all NaNs.
        hasNaNs = any(isnan(T(:)));
        hasInfs = any(isinf(T(:)));
        if ( hasNaNs || hasInfs )
            R1 = tril(NaN(size(T)));
            return;
        end
        
        % (4) We compute chol(T,'lower'). If this errors, we compute 
        % chol(T + delta*eye(size(T)),'lower') where delta is a small 
        % number like eps(class(T)). If this also errors out, we keep 
        % increasing delta until we get no errors.
            epsilon = eps(class(T));
            delta = epsilon;
            I = eye(size(T));
            found = false;
            while( found == false )
                [R1,p] = chol(T + delta*I,'lower');
                if ( p == 0 )
                    % We are done.
                    found = true;
                else
                    % Increase delta.
                    %delta = delta + epsilon;
                    delta = 2*delta;
                end
            end % end of while.
            
        end % end of singularLowerChol.        
        
        %=== helper function createCovariance
        function cmat = createCovariance(name,dimension,varargin)
%createCovariance - Create a specified covariance matrix.
%   CMAT = createCovariance(NAME,DIMENSION,varargin) takes a string NAME,
%   and a positive integer DIMENSION and optional arguments and creates a
%   covariance matrix corresponding to NAME of the specified DIMENSION.
%   Permissible values of NAME are as follows:
%
%       TYPE_FULL = 'Full';
%       TYPE_FULLCHOLESKY = 'FullCholesky';
%       TYPE_DIAGONAL = 'Diagonal';
%       TYPE_ISOTROPIC = 'Isotropic';
%       TYPE_COMPSYMM = 'CompSymm';
%       TYPE_FIXEDWEIGHTS = 'FixedWeights';
%       TYPE_PATTERNED = 'Patterned';
%
%   We don't allow BlockedCovariance to be constructed here. For that, use
%   the BlockedCovariance constructor.

            % (1) Make sure we have at least 2 input arguments.
            narginchk(2,Inf);
            
            % (2) Ensure that name is sensible.
            name = internal.stats.getParamVal(name,classreg.regr.lmeutils.covmats.CovarianceMatrix.AllowedCovarianceTypes,'NAME');                        
            
            % (3) Construct the appropriate object.
            switch lower(name)
                case lower(classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FULL)                    
                    cmat = classreg.regr.lmeutils.covmats.FullCovariance(dimension,varargin{:});    
                    
                case lower(classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FULLCHOLESKY)
                    cmat = classreg.regr.lmeutils.covmats.FullCholeskyCovariance(dimension,varargin{:});
                    
                case lower(classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_DIAGONAL)
                    cmat = classreg.regr.lmeutils.covmats.DiagonalCovariance(dimension,varargin{:});
                    
                case lower(classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_ISOTROPIC)
                    cmat = classreg.regr.lmeutils.covmats.IsotropicCovariance(dimension,varargin{:});
                    
                case lower(classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_COMPSYMM)
                    cmat = classreg.regr.lmeutils.covmats.CompoundSymmetryCovariance(dimension,varargin{:});                                    
                    
                case lower(classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FIXEDWEIGHTS)
                    cmat = classreg.regr.lmeutils.covmats.FixedWeightsCovariance(dimension,varargin{:});
                    
                case lower(classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_PATTERNED)                    
                    cmat = classreg.regr.lmeutils.covmats.PatternedCovariance(dimension,varargin{:});
                    
                otherwise
                    % <entry key="BadCovarianceName">Covariance {0} is not currently supported.</entry>
                    error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadCovarianceName',name));                    
            end
            
        end % end of createCovariance.
        
    end
    
end

%=== helper method for verifying assertions
function assertThat(condition,msgID,varargin)
%   assertThat(condition,msgID,varargin) takes a variable condition that is
%   either true or false, a message catalog ID msgID and optional arguments
%   required to create a message object from msgID. If condition is false,
%   an error message represented by msgID is thrown.

    if ~condition        
        % (1) Create a message object from msgID and varargin.
        try
            msg = message(msgID,varargin{:});
        catch
            % <entry key="BadMsgID">Invalid message ID: {0}.</entry>      
            error(message('stats:LinearMixedModel:BadMsgID',msgID));
        end        
        % (2) Create and throw an MException.
        ME = MException( msg.Identifier, getString(msg) );
        throwAsCaller(ME);        
    end

end % end of assertThat.

