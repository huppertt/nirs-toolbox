classdef BlockedCovariance < classreg.regr.lmeutils.covmats.CovarianceMatrix
%  This feature is intended for internal use only and is subject to change
%  at any time without warning. 

%BlockedCovariance - Class to represent a 'Blocked' covariance matrix.
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
%   When PSI is of Type 'Blocked', D has a block diagonal structure.
%   Each block of D is optionally itself block diagonal. For example,
%   consider the following D:
%
%     D = [ D_1   |
%              D_1|
%           ------|----------
%                 |D_2      |
%                 |   D_2   |
%                 |      D_2|
%                 |---------|------|
%                           |D_3   |
%                           |   D_3|
%                           |------| ]
%
%   Here D is block diagonal with 3 blocks. The first block is itself block
%   diagonal and is made up of q_1 by q_1 matrices D_1 with 2 repetitions. 
%   The second block is also block diagonal and is made up of q_2 by q_2 
%   matrices D_2 with 3 repetitions. The third block is also block diagonal
%   and is made up of q_3 by q_3 matrices D_3 with 2 repetitions. If
%   theta_1, theta_2 and theta_3 are the unconstrained parameter vectors
%   for D_1, D_2 and D_3 respectively then 
%
%            theta = [theta_1;theta_2;theta_3] 
%
%   is the unconstrained parameter vector for D. The Cholesky factor of D 
%   is given by:
%
%     L = [ L_1   |
%              L_1|
%           ------|----------
%                 |L_2      |
%                 |   L_2   |
%                 |      L_2|
%                 |---------|------|
%                           |L_3   |
%                           |   L_3|
%                           |------| ]
%
%   where L_1, L_2 and L_3 are the Cholesky factors of D_1, D_2 and D_3
%   respectively.
%
%   (1) The unconstrained parameterization of PSI in this case is given by 
%       the unconstrained parameters theta = [theta_1;theta_2;theta_3] and 
%       log(sigma).
%
%   (2) If eta_1, eta_2 and eta_3 are the natural parameter vectors for 
%       sigma^2 * D_1, sigma^2 * D_2 and sigma^2 * D_3 then eta =
%       [eta_1;eta_2;eta_3] will be the natural parameter vector for
%       sigma^2 * D. The natural parameterization of PSI is then given by
%       the natural parameters eta and log(sigma).
%
%   (3) If heta_1, heta_2 and heta_3 are the canonical parameter vectors
%       for sigma^2 * D_1, sigma^2 * D_2 and sigma^2 * D_3 then heta =
%       [heta1;heta2;heta3] will be the canonical parameter vector for
%       sigma^2 * D. The canonical parameterization of PSI is then given by
%       the canonical parameters heta and sigma.
%
%   BlockedCovariance methods:
%       BlockedCovariance - Create a 'Blocked' covariance matrix.
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
%   BlockedCovariance properties:
%       Name - Name of this covariance matrix
%       Size - Size of this covariance matrix
%       Type - The type of this covariance matrix (e.g., "Diagonal").
%       CovariancePattern - Pattern of non-zeros in this covariance matrix
%       WhichParameterization - Parameterization that is currently in use.
%       NumBlocks - Number of blocks in this covariance matrix.
%       Matrices - Objects for the smaller matrices in BlockedCovariance.
%       NumReps - Number of repetitions of smaller matrix in each block.
%       SizeVec - Size of smaller matrices that make up BlockedCovariance.
%       NumParametersExcludingSigmaVec - Parameters in matrices that form BlockCovariance.
%
%   See also classreg.regr.lmeutils.covmats.CovarianceMatrix
    
%   Copyright 2012-2013 The MathWorks, Inc.

    properties (GetAccess=public, SetAccess=protected)                
%NumBlocks - The number of blocks in this BlockedCovariance. Each block is
%   optionally itself block diagonal.
        NumBlocks

%Matrices - A cell vector of size NumBlocks by 1 containing objects
%   representing the smaller matrices that make up this BlockedCovariance.
%   Block i of BlockedCovariance contains NumReps(i) repetitions of the
%   matrix represented by Matrices{i}.
        Matrices       
        
%NumReps - A vector of size NumBlocks by 1 indicating the number of
%   repetitions of the smaller matrices forming each block.
        NumReps       
        
%SizeVec - A vector containing the size of smaller matrices that make up
%   this BlockedCovariance. SizeVec(i) is equal to Matrices{i}.Size.
        SizeVec              

%NumParametersExcludingSigmaVec - The number of unconstrained parameters 
%   in the smaller matrices that make up this BlockedCovariance. 
%   NumParametersExcludingSigmaVec(i) is equal to Matrices{i}.NumParametersExcludingSigma.
        NumParametersExcludingSigmaVec                 
    end

    properties (Access=private)                
        % Linear index into non-zero elements based on CovariancePattern
        D_nzind
        
        % Linear index into non-zero elements of tril(CovariancePattern)
        L_nzind
        
        % Lower triangular part of CovariancePattern
        Lpat        
    end
    
    methods (Access=public)
        
        function this = BlockedCovariance(mats,reps,varargin)
%

%BlockedCovariance - Construct a BlockedCovariance object.
%   OBJ = BlockedCovariance(MATS,REPS) constructs an object OBJ of type
%   BlockedCovariance from a cell array of CovarianceMatrix objects MATS 
%   and a vector of integers REPS. MATS and REPS must have the same length.
%   The BlockedCovariance OBJ will have length(MATS) blocks. Block i will 
%   be made up of REPS(i) repetitions of matrix represented by MATS{i} 
%   along the block diagonal.           
%                        
%   OBJ = BlockedCovarianceMatrix(...,'Name',NAME) also specifies the name to be 
%   associated with this covariance matrix.

           this = this@classreg.regr.lmeutils.covmats.CovarianceMatrix();

            % Process input args
            narginchk(2,Inf);                                                            
            dfltName = [];
            parnames = {  'Name'  };
               dflts = { dfltName};                            
            name = internal.stats.parseArgs(parnames,dflts,varargin{:});                                                
            
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
               name = 'Blocked Covariance Matrix';
            end
  
            % Validate mats.
            if ~iscell(mats)
                % <entry key="BadMatrices">MATS must be a cell array of CovarianceMatrix objects.</entry>                     
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadMatrices'));                
            end
            
            % Make mats into a column vector
            if size(mats,1) == 1
                mats = mats';
            end
            
            % Make sure reps is a vector of positive integers
            if ~internal.stats.isIntegerVals(reps,1)
                % <entry key="BadReps">REPS must be a vector of positive integers.</entry>
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadReps'));
            end
            
            % Make reps into a column vector
            if size(reps,1) == 1
                reps = reps';
            end
            
            % Make sure mats and reps have the same length
            if length(mats) ~= length(reps)
                % <entry key="InvalidMatricesReps">MATS and REPS must be vectors of the same length.</entry>
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:InvalidMatricesReps'));                
            end
            
            % Ensure elements of mats are subclasses of CovarianceMatrix
            numblocks = length(mats);
            for i = 1:numblocks               
               if ~isa(mats{i},'classreg.regr.lmeutils.covmats.CovarianceMatrix')                   
                   % <entry key="BadMatrices">MATS must be a cell array of CovarianceMatrix objects.</entry>                     
                   error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:BadMatrices'));                   
               end
            end
            
            % Set number of blocks.            
            this.NumBlocks = numblocks;
            
            % Set component matrices.
            this.Matrices = mats;
            
            % Set repetitions of component matrices.                        
            this.NumReps = reps;
            
            % Set name for this covariance matrix.
            this.Name = name;
                                    
            % Set Size and SizeVec
            sizevec = zeros(numblocks,1);
            for i = 1:numblocks
                sizevec(i) = mats{i}.Size;
            end
            this.Size = sum(sizevec.*reps);
            this.SizeVec = sizevec;
                          
            % Set NumParameters and NumParametersvec
            numparameters = 0;
            numparametersvec = zeros(numblocks,1);
            for i = 1:numblocks
               numparametersvec(i) = mats{i}.NumParametersExcludingSigma;
               numparameters = numparameters + numparametersvec(i); 
            end
            this.NumParametersExcludingSigma = numparameters;
            this.NumParametersExcludingSigmaVec = numparametersvec;
            
            % VariableNames does not make sense for a BlockedCovariance.
            % So set it to [].            
            this.VariableNames = [];                        
            
            % Use implementations of abstract methods inherited from
            % CovarianceMatrix to set the Type, CovariancePattern,
            % WhichParameterization, UnconstrainedParameters,
            % NumParametersExcludingSigma, CanonicalParameterNames and
            % Sigma.            
            this.Type = defineType(this);      
            assert( internal.stats.isString(this.Type,false) == true );
            
            this.CovariancePattern = defineCovariancePattern(this);
            dimension = this.Size;
            assert( all(size(this.CovariancePattern) == [dimension,dimension]) );            
            
            this.WhichParameterization = classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_UNCONSTRAINED;
            theta = initializeTheta(this);
            this.UnconstrainedParameters = theta;
            
            this.CanonicalParameterNames = defineCanonicalParameterNames(this);
            
            this.Sigma = 1;                                                
            
            
            % Set private properties: D_nzind, Lpat and L_nzind            
            this.D_nzind = find(this.CovariancePattern);
            this.Lpat = tril(this.CovariancePattern);
            this.L_nzind = find(this.Lpat);     
            
        end % end of BlockedCovariance.
        
        function verifyTransformations(this)
%verifyTransformations - Check the parameter transformations between 
%   unconstrained, natural and canonical for accuracy.
          
            % (1) theta to D and D to theta
            theta = randn(this.NumParametersExcludingSigma,1);
            D = theta2D(this,theta);            
            theta_recon = D2theta(this,D);            
            err(1) = max(abs(theta(:) - theta_recon(:)));            
           
            % (2) 
            % (A) L via chol vs. getLowerTriangularCholeskyFactor
            theta = randn(this.NumParametersExcludingSigma,1);
            D = theta2D(this,theta);            
            L = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(D);            
            this = setUnconstrainedParameters(this,theta);
            L1 = getLowerTriangularCholeskyFactor(this);            
            err(2) = max(abs(L(:) - L1(:)));

            % (3)
            % (B) theta to natural and then getLowerTriangularCholeskyFactor
            natural = convertUnconstrainedToNatural(this,theta);
            this = setNaturalParameters(this,natural);
            L2 = getLowerTriangularCholeskyFactor(this);            
            err(3) = max(abs(L(:) - L2(:)));
            
            % (3) 
            % (C) canonical parameters from unconstrained vs. natural
            this = setUnconstrainedParameters(this,theta);
            canonical1 = getCanonicalParameters(this);
            this = setNaturalParameters(this,natural);
            canonical2 = getCanonicalParameters(this);
            err(4) = max(abs(canonical1(:) - canonical2(:)));
                                    
            if max(err(:)) <= sqrt(eps)
                disp(getString(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:String_TransformsOK')));                
            else
                error(message('stats:classreg:regr:lmeutils:covmats:CovarianceMatrix:String_TransformsNotOK'));                
            end            
            
        end % end of verifyTransformations.       
        
        function this = setUnconstrainedParameters(this,theta)
%setUnconstrainedParameters - Sets the unconstrained parameterization.
%   OBJ = setUnconstrainedParameters(OBJ,THETA) accepts an object OBJ of 
%   type CovarianceMatrix, an unconstrained parameter vector THETA and 
%   returns the object OBJ with the specified unconstrained parameters. 
%   The length of THETA must match the number of unconstrained parameters
%   (excluding sigma) in OBJ.

            this = setUnconstrainedParameters@classreg.regr.lmeutils.covmats.CovarianceMatrix(this,theta);
            
            numblocks = this.NumBlocks;
            offset = 0;
            for i = 1:numblocks
                startidx = offset + 1;
                endidx = offset + this.NumParametersExcludingSigmaVec(i);
                thetai = theta(startidx:endidx);
                offset = endidx;
                this.Matrices{i} = this.Matrices{i}.setUnconstrainedParameters(thetai);
                this.Matrices{i}.WhichParameterization = ...
                    classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_UNCONSTRAINED;
            end                        

        end
        
        function this = setNaturalParameters(this,eta)
%setNaturalParameters - Sets the Natural parameterization.
%   OBJ = setNaturalParameters(OBJ,ETA) accepts an object OBJ of type
%   CovarianceMatrix, a Natural parameter vector ETA and returns the object
%   OBJ with the specified Natural parameters. The length of ETA must match
%   the number of unconstrained parameters (excluding sigma) in OBJ.

            this = setNaturalParameters@classreg.regr.lmeutils.covmats.CovarianceMatrix(this,eta);
            
            numblocks = this.NumBlocks;
            offset = 0;
            for i = 1:numblocks
                startidx = offset + 1;
                endidx = offset + this.NumParametersExcludingSigmaVec(i);
                etai = eta(startidx:endidx);
                offset = endidx;
                this.Matrices{i} = this.Matrices{i}.setNaturalParameters(etai);
                this.Matrices{i}.WhichParameterization = ...
                    classreg.regr.lmeutils.covmats.CovarianceMatrix.PARAMETERIZATION_NATURAL;
            end

        end
        
         function this = setSigma(this,sigma)
%setSigma - Set sigma for this covariance matrix.
%   OBJ = setSigma(OBJ,SIGMA) accepts an object of type CovarianceMatrix, 
%   a real positive scalar SIGMA, sets the Sigma property of OBJ and 
%   returns the modified object.

            this = setSigma@classreg.regr.lmeutils.covmats.CovarianceMatrix(this,sigma);
            
            numblocks = this.NumBlocks;
            for i = 1:numblocks
                this.Matrices{i} = this.Matrices{i}.setSigma(sigma);
            end
            
         end
        
    end
    
    methods (Access=public, Hidden=true)
        
        % A subclass must know how to initialize theta and how to convert
        % to L. For verifying transforms, it is also useful for subclasses
        % to define theta2D and D2theta methods.
        function theta = initializeTheta(this)
            
            theta = zeros(this.NumParametersExcludingSigma,1);
            offset = 0;
            for i = 1:this.NumBlocks
                startidx = offset + 1;
                endidx = offset + this.NumParametersExcludingSigmaVec(i);
                theta(startidx:endidx) = this.Matrices{i}.initializeTheta();
                offset = endidx;
            end
                        
        end % end of initializeTheta.        
        
        function L = theta2L(this,theta)            
                       
            % Get number of blocks.
            numblocks = this.NumBlocks;
            
            % Get vectors l_1,l_2,...,l_r from L_1,L_2,..,L_r.
            % Replicate l_i, NumReps(i) times and store in lvec{i}.
            lvec = cell(numblocks,1);
            offset = 0;
            for i = 1:numblocks
                startidx = offset + 1;
                endidx = offset + this.NumParametersExcludingSigmaVec(i);                                
                thetai = theta(startidx:endidx);
                offset = endidx;
                Li = this.Matrices{i}.theta2L(thetai);
                % Store non-zero elements of lower triangular Li into
                % vector li.
                li = Li( tril(this.Matrices{i}.CovariancePattern) );
                %lvec{i} = kron( ones(this.NumReps(i),1), li );
                lvec{i} = repmat(li,this.NumReps(i),1);
            end
            
            % L contains 1's in places where there should be a non-zero
            % and L_nzind contains a linear index into nonzero elements
            % in L (columnwise order).
            L = double(this.Lpat);
            L(this.L_nzind) = cell2mat(lvec(:));
                                                
        end % end of theta2L.
        
        function D = theta2D(this,theta)
                                   
            % Get number of blocks.
            numblocks = this.NumBlocks;
            
            % Get vectors d_1,d_2,...,d_r from D_1,D_2,..,D_r.
            % Replicate d_i, NumReps(i) times and store in dvec{i}.
            dvec = cell(numblocks,1);
            offset = 0;
            for i = 1:numblocks
                startidx = offset + 1;
                endidx = offset + this.NumParametersExcludingSigmaVec(i);                                
                thetai = theta(startidx:endidx);
                offset = endidx;
                Di = this.Matrices{i}.theta2D(thetai);
                % Store non-zero elements of Di in a vector di
                di = Di(this.Matrices{i}.CovariancePattern);
                dvec{i} = kron( ones(this.NumReps(i),1), di );
            end
            
            % D contains 1's in places where there should be a non-zero
            % and D_nzind contains a linear index into nonzero elements
            % in D (columnwise order).
            D = double(this.CovariancePattern);
            D(this.D_nzind) = cell2mat(dvec(:));
            
        end % end of theta2D.
        
        function theta = D2theta(this,D)
            
            theta = zeros(this.NumParametersExcludingSigma,1);
            offset = 0;
            for i = 1:this.NumBlocks
                startidx = offset + 1;
                endidx = offset + this.NumParametersExcludingSigmaVec(i);                                
                % Look at the 1st matrix in the ith block.
                Dsub = getSubMatrix(this,D,i,1);
                theta(startidx:endidx) = this.Matrices{i}.D2theta(Dsub);
                offset = endidx;
            end
                        
        end % end of D2theta.
               
        % A subclass must know how to convert between unconstrained and
        % natural parameterizations.
        function natural = convertUnconstrainedToNatural(this,unconstrained)
            
            assert(length(unconstrained) == this.NumParametersExcludingSigma);
            
            natural = zeros(this.NumParametersExcludingSigma,1);
            offset = 0;
            for i = 1:this.NumBlocks
                startidx = offset + 1;
                endidx = offset + this.NumParametersExcludingSigmaVec(i);                                
                natural(startidx:endidx) = ...
                    this.Matrices{i}.convertUnconstrainedToNatural(unconstrained(startidx:endidx));                
                offset = endidx;
            end
            
            
        end % end of convertUnconstrainedToNatural.
                
        function unconstrained = convertNaturalToUnconstrained(this,natural)
            
            assert(length(natural) == this.NumParametersExcludingSigma);
            
            unconstrained = zeros(this.NumParametersExcludingSigma,1);
            offset = 0;
            for i = 1:this.NumBlocks
                startidx = offset + 1;
                endidx = offset + this.NumParametersExcludingSigmaVec(i);                                
                unconstrained(startidx:endidx) = ...
                    this.Matrices{i}.convertNaturalToUnconstrained(natural(startidx:endidx));                
                offset = endidx;
            end
            
        end % end of convertNaturalToUnconstrained.
        
        % A subclass must know how to convert between natural and canonical
        % parameterizations.
        function canonical = convertNaturalToCanonical(this,natural)
           
            assert(length(natural) == this.NumParametersExcludingSigma);
            
            canonical = zeros(this.NumParametersExcludingSigma,1);
            offset = 0;
            for i = 1:this.NumBlocks
                startidx = offset + 1;
                endidx = offset + this.NumParametersExcludingSigmaVec(i);                                
                canonical(startidx:endidx) = ...
                    this.Matrices{i}.convertNaturalToCanonical(natural(startidx:endidx));                
                offset = endidx;
            end
            
        end % end of convertNaturalToCanonical.
        
        function natural = convertCanonicalToNatural(this,canonical)
            
            assert(length(canonical) == this.NumParametersExcludingSigma);
            
            natural = zeros(this.NumParametersExcludingSigma,1);
            offset = 0;
            for i = 1:this.NumBlocks
                startidx = offset + 1;
                endidx = offset + this.NumParametersExcludingSigmaVec(i);                                
                natural(startidx:endidx) = ...
                    this.Matrices{i}.convertCanonicalToNatural(canonical(startidx:endidx));                
                offset = endidx;
            end
            
        end % end of convertCanonicalToNatural.
        
        % A subclass must define canonical parameter names.
        function names = defineCanonicalParameterNames(this)
           
            % Loop over blocks. Get the canonical parameter names for each
            % child matrix and save it into an element of the cell array
            % names.
            numblocks = this.NumBlocks;
            names = cell(numblocks,1);
            for i = 1:numblocks
                names{i} = this.Matrices{i}.defineCanonicalParameterNames;
            end
            
        end % end of defineCanonicalParameterNames.
        
        % A subclass must define its type (e.g., 'Full','Diagonal' etc.)
        function type = defineType(this) %#ok<MANU>
        % A 'Blocked' covariance matrix has type Blocked. 
            type = classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_BLOCKED;
        end % end of defineType.
        
        % A subclass must define its covariance pattern.
        function covariancepattern = defineCovariancePattern(this)

            % Loop over blocks. Within each block loop over reps of the
            % covariance matrix in the block. Calculate offset from top
            % left corner into the covariance pattern and set it using the
            % covariance pattern of the current block.
            numblocks = this.NumBlocks;
            reps = this.NumReps;
            localq = this.SizeVec;            
            % Form the product localq * reps (element by element)
            qreps = localq .* reps;            
            covariancepattern = logical(sparse( sum(qreps), sum(qreps) ));
            for r = 1:numblocks              
                for k = 1:reps(r)
                    offset = sum(qreps(1:(r-1))) + (k-1)*localq(r);
                    idx = offset + 1 : offset + localq(r);
                    covariancepattern(idx,idx) = this.Matrices{r}.CovariancePattern; 
                end
            end           
            
        end % end of defineCovariancePattern.
        
        
    end
    
    methods (Access=private)
       
        function Dsub = getSubMatrix(this,D,block,rep)
%getSubMatrix - Get submatrix of big matrix D.
%   Dsub = getSubMatrix(this,D,block,rep) gets the submatrix of D located 
%   at block number block and rep number rep. No input validation.
    
            % Block 1 has size qvec(1)*NumReps(1), Block 2 has size
            % qvec(2)*NumReps(2) etc.
            qreps = this.SizeVec .* this.NumReps;

            % Skip (block-1) blocks and (rep-1) reps (size SizeVec(block))
            offset = sum(qreps(1:block-1)) + (rep-1)*this.SizeVec(block);
            
            % row and column indices that we want
            idx = offset + 1 : offset + this.SizeVec(block);
            
            % Get Dsub
            Dsub = D(idx,idx);
        end % end of getSubMatrix.
                
    end
    
end


