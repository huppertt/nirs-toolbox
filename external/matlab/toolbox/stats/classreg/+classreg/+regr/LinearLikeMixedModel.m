classdef (Abstract) LinearLikeMixedModel < classreg.regr.ParametricRegression
%  This feature is intended for internal use only and is subject to change
%  at any time without warning.

%LinearLikeMixedModel - An abstract class that defines common properties 
%   and methods shared by GeneralizedLinearMixedModel and LinearMixedModel
%   which are the concrete subclasses of this class.
    
% Shared internal constants.
properties (Constant=true,Hidden=true)
%AllowedDummyVarCodings - Dummy variable codings for categorical variables 
%   that we currently support.       
        AllowedDummyVarCodings = {'reference','referencelast','effects',...
            'full','difference','backwardDifference'};
    
%AllowedCovariancePatterns - Covariance patterns for the random effects 
%   vector that we currently support.
        AllowedCovariancePatterns = {classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FULL,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FULLCHOLESKY,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_DIAGONAL,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_ISOTROPIC,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_COMPSYMM,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FIXEDWEIGHTS,...
            classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_PATTERNED};    
end

% Shared properties storing design matrix info and a standard form model.
properties (Access={?classreg.regr.LinearLikeMixedModel})
%   The following variables only include observations used in the fit 
%   corresponding to ObservationInfo.Subset.    

%y - Response vector for the model.
        y
        
%FixedInfo - Structure storing fixed effects info. See extractFixedInfo. 
%   FixedInfo has the following fields:
%
%             X = Fixed Effects design matrix.
%
%     XColNames = Cell array containing names of columns of X.
%
%   XCols2Terms = Integer array storing which term a particular column
%                 of X belongs to. XCols2Terms(i) = j means that the
%                 ith column belongs to "term" j.
        FixedInfo    
       
%RandomInfo - Structure storing random effects info. See extractRandomInfo.
%   RandomInfo has the following fields:
%
%             Z = R-by-1 cell array of random effects design matrices.
%
%             q = R-by-1 vector with q(i) = size(Z{i},2).
%
%     ZColNames = R-by-1 cell array. ZColNames{i} is a cell array that 
%                 contains the column names of Z{i}.
%
%    ZsColNames = Table containing the name of the columns of sparse
%                 matrix Zs created in fitter. This same Zs is also stored
%                 in object slme of type StandardLinearMixedModel. Row i of
%                 ZsColNames names column i of Zs. ZsColNames is a table
%                 with 3 variables:
%
%                       Group   Level   Name
%
%                 Group = group associated with this random effect
%                 Level = level associated with this random effect
%                 Name  = predictor name associated with this random effect
        RandomInfo    
        
%GroupingInfo - Structure storing grouping variable info. See extractGroupingInfo.
%   GroupingInfo has the following fields:
%
%             R = number of random effects grouping variables.
%
%             G = R-by-1 cell array of grouping variables supplied by
%                 the user. G{i} is a grouping variable for Z{i}.     
%
%        GNames = R-by-1 cell array containing the names of grouping 
%                 variables. GNames{i} is the name of the grouping 
%                 variable G{i}. 
%
%           Gid = R-by-1 cell array of grouping variables after passing
%                 each G{i} through grp2idx:
%                 [Gid{i},GidLevelNames{i}] = grp2idx(G{i});
%
% GidLevelNames = A R-by-1 cell array such that GidLevelNames{i}{j} is 
%                 mapped to "ID" j for Gid{i}.    
%
%           lev = R-by-1 array containing the number of levels in each
%                 grouping variable. lev(i) is the length of 
%                 GidLevelName{i} which is also the number of levels in
%                 grouping variable Gid{i}.    
        GroupingInfo
    
%slme - A StandardLinearMixedModel object for a LinearMixedModel and a 
%   StandardGeneralizedLinearMixedModel object for a
%   GeneralizedLinearMixedModel.
        slme
end

% Shared utility methods.
methods (Access={?classreg.regr.LinearLikeMixedModel})
    
    function y = extractResponse(model,ds)
%extractResponse - Extract the response variable from a table ds.
%   y = extractResponse(model,ds) accepts a LinearLikeMixedModel object 
%   model and a table ds and extracts the response vector.

        % (1) Extract response var from ds.
        responseName = model.Formula.ResponseName;
        y = ds.(responseName);        
        
        % (2) Ensure that y is a numeric real vector.
        assertThat(isnumeric(y) & isreal(y) & isvector(y),'stats:LinearMixedModel:BadResponseVar',responseName);
        
    end % extractResponse.

    function FixedInfo = extractFixedInfo(model,ds)
%extractFixedInfo - Extracts info about fixed effects from a table ds.    
%   FixedInfo = extractFixedInfo(model,ds) accepts a LinearLikeMixedModel 
%   object model and a table ds and extracts the structure FixedInfo.
%
%   FixedInfo has the following fields:
%
%             X = Fixed Effects design matrix.
%
%     XColNames = Names of columns of X.
%
%   XCols2Terms = Integer array storing which term a particular column
%                 of X belongs to. XCols2Terms(i) = j means that the
%                 ith column belongs to "term" j. 
        
        [FixedInfo.X,FixedInfo.XColNames,FixedInfo.XCols2Terms] = ...
            createDesignMatrix(model,ds,model.Formula.FELinearFormula);
        
    end % end of extractFixedInfo.

    function RandomInfo = extractRandomInfo(model,ds)
%extractRandomInfo - Extracts info about random effects from a table ds.
%   RandomInfo = extractRandomInfo(model,ds) accepts a LinearLikeMixedModel 
%   object model and a table ds and extracts the structure RandomInfo. 
%
%   RandomInfo has the following fields:
%
%             Z = R-by-1 cell array of random effects design matrices.
%
%             q = R-by-1 vector with q(i) = size(Z{i},2).
%
%     ZColNames = R-by-1 cell array. ZColNames{i} contains the column
%                 names of Z{i}.

        % (1) Call designFromFormula.
        R = length(model.Formula.RELinearFormula);
        RandomInfo.Z = cell(R,1);
        RandomInfo.ZColNames = cell(R,1);
        RandomInfo.q = zeros(R,1);
        for i = 1:R             
            [RandomInfo.Z{i},RandomInfo.ZColNames{i},~] = ...
                createDesignMatrix(model,ds,model.Formula.RELinearFormula{i});                
            RandomInfo.q(i) = size(RandomInfo.Z{i},2);
        end                        
        
    end % end of extractRandomInfo.

    function GroupingInfo = extractGroupingInfo(model,ds)
%extractGroupingInfo - Extracts grouping info from a table ds.
%   GroupingInfo = extractGroupingInfo(model,ds) accepts a
%   LinearLikeMixedModel object model and a table ds and extracts a
%   structure GroupingInfo.
%
%   GroupingInfo has the following fields:
%
%             R = number of random effects grouping variables.
%
%             G = R-by-1 cell array of grouping variables supplied by
%                 the user. G{i} is a grouping variable for Z{i}.     
%
%        GNames = R-by-1 cell array containing the names of grouping 
%                 variables. GNames{i} is the name of the grouping 
%                 variable G{i}. 
%
%           Gid = R-by-1 cell array of grouping variables after passing
%                 each G{i} through grp2idx:
%                 [Gid{i},GidLevelNames{i}] = grp2idx(G{i});
%
% GidLevelNames = A R-by-1 cell array such that GidLevelNames{i}{j} is 
%                 mapped to "ID" j for Gid{i}.    
%
%           lev = R-by-1 array containing the number of levels in each
%                 grouping variable. lev(i) is the number of levels in
%                 grouping variable Gid{i}.    

        % (1) First make G and GNames.
        R = length(model.Formula.RELinearFormula);
        G = cell(R,1);
        GNames = cell(R,1);
        for i = 1:R           
            % Which variables make up this interaction?
            interactionVars = model.Formula.GroupingVariableNames{i};
            
            % Combine interactionVars into a single grouping variable.
            [G{i},GNames{i}] = ...
                model.makeInteractionVar(ds,interactionVars);                           
        end

        % (2) Pass each G{i} through grp2idx to get Gid and GidLevelNames.
        Gid = cell(R,1);
        GidLevelNames = cell(R,1);
        lev = zeros(R,1);
        for i = 1:R
            [Gid{i},GidLevelNames{i}] = grp2idx(G{i});
            lev(i) = length(GidLevelNames{i});
            %lev(i) = length(unique(Gid{i}));
        end        
        
        % (3) Store everything in GroupingInfo.
        GroupingInfo.R = R;
        GroupingInfo.G = G;
        GroupingInfo.GNames = GNames;
        GroupingInfo.Gid = Gid;
        GroupingInfo.GidLevelNames = GidLevelNames;
        GroupingInfo.lev = lev;
        
    end % end of extractGroupingInfo.

    function Psi = makeCovarianceMatrix(model)
%makeCovarianceMatrix - Makes a blocked covariance object. 
%   Psi = makeCovarianceMatrix(model) accepts a LinearLikeMixedModel object 
%   model and returns Psi, an object of type BlockedCovariance. Psi
%   represents the covariance matrix of the model in standard form.
%
%   covariancepattern{i} is a string (such as 'Full' or 'Diagonal')
%   representing the structure of the covariance matrix for grouping
%   variable i. The length of covariancepattern must be equal to
%   model.GroupingInfo.R.

        % (1) Make component submatrices for the BlockedCovariance object.
        R = model.GroupingInfo.R;        
        covariancepattern = model.CovariancePattern;
        assert(R == length(covariancepattern));        
        mat = cell(R,1);
        for i = 1:R     

            % (2) covariancepattern{i} may be a string or a logical matrix.
            % If it is a logical matrix then we make a PatternCovariance
            % object with the covariance pattern in covariancepattern{i}.
            if islogical(covariancepattern{i})
                mat{i} = classreg.regr.lmeutils.covmats.CovarianceMatrix.createCovariance(classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_PATTERNED,...
                    model.RandomInfo.q(i),'Name',model.GroupingInfo.GNames{i},...
                    'VariableNames',model.RandomInfo.ZColNames{i},'CovariancePattern',covariancepattern{i});
            else
                mat{i} = classreg.regr.lmeutils.covmats.CovarianceMatrix.createCovariance(covariancepattern{i},...
                    model.RandomInfo.q(i),'Name',model.GroupingInfo.GNames{i},...
                    'VariableNames',model.RandomInfo.ZColNames{i});
            end
            
        end
        
        % (3) Make the BlockedCovariance object.
        if ( R == 0 )                            
            Psi = classreg.regr.lmeutils.covmats.BlockedCovariance({classreg.regr.lmeutils.covmats.FullCovariance(0)},1);
        else
            Psi = classreg.regr.lmeutils.covmats.BlockedCovariance(mat,model.GroupingInfo.lev);
        end
        
    end % end of makeCovarianceMatrix.

    function ZsColNames = makeSparseZNames(model)
%makeSparseZNames - Make names for the columns of big sparse Z matrix.        
%   ZsColNames = makeSparseZNames(model) takes an object model of type
%   LinearLikeMixedModel and makes ZsColNames, a table array such that each 
%   column of Zs returned by makeSparseZ is named by the corresponding row 
%   in ZsColNames. Every level of every grouping variable gets a random 
%   effects vector of the right size.
%
%   NOTE:
%
%   Z = R-by-1 cell array. Z{i} = N-by-q(i) design matrix for grouping variable i.
%   q = R-by-1 vector with q(i) = length of random effect vector for ith grouping variable.
% lev = R-by-1 array. lev(i) = number of levels in the ith grouping variable.
% Gid = R-by-1 cell array. Gid{i} = N-by-1 vector containing integers from 1 to lev(i)
%       indicating the membership of each observation into one of the levels from 1
%       to lev(i) of grouping variable i.
%
%%
% Suppose there are $r$ grouping variables. For the $i$ th grouping variable, we
% have $lev(i)$ levels. In each level of the $i$ th grouping variable, there is
% a random effect vector of length $q(i)$ which can be considered a draw from
% $N(0,\sigma^2 D_i)$. The co-variance $D_i$ is parameterized by a vector
% $\theta_i$ of length $lenTheta(i)$. The input vector $\theta$ is a concatenation
% of such vectors $\theta_i$. Thus:
%
% $\theta = [\theta_1; \theta_2; ...; \theta_r]$ and
%
% $length(\theta_i) = lenTheta(i)$
%
% If we concatenate all random effect vectors across all grouping variables
% into one big vector, we get the following vector:
%
% $b =
% [b(1)_1,b(1)_2,..,b(1)_{lev(1)},b(2)_1,b(2)_2,..,b(2)_{lev(2)}...,b(r)_1,b(r)_2,..,b(r)_{lev(r)}]^T$
%
% $b(i)_j = j$ th level random effect vector for $i$ th grouping variable
% of length $q(i) \times 1$
%
% The length of $b$ will be $\sum_{i = 1}^r q(i) lev(i)$.
%
% To understand, how Zs is constructed from Z, let us take a simple
% example. Suppose:
%
% Z{1} =
%
%     1     2
%     3     4
%     5     6
%     7     8
%     9    10
%
% Z{2} =
%
%    11    12    13
%    14    15    16
%    17    18    19
%    20    21    22
%    23    24    25
%
% lev(1) =
%
%     2
%
% lev(2) =
%
%     3
%
% Gid{1} =
%
%     1
%     2
%     1
%     1
%     2
%
% Gid{2} =
%
%     1
%     2
%     3
%     3
%     2
%
% The columns in matrix Zs correspond to the concatenated $b$ vector:
% $b =
% [b(1)_1,b(1)_2,..,b(1)_{lev(1)},b(2)_1,b(2)_2,..,b(2)_{lev(2)}...,b(r)_1,b(r)_2,..,b(r)_{lev(r)}]^T$
%
% All rows of Z{1} that satisfy Gid{1} == 1 go into columns 1:q(1) of Zs.
% If idx = Gid{1} == 1 then
% Zs(idx,1:q(1)) = Z{1}(idx,:)
%
% All rows of Z{r} that satisfy Gid{r} == k go into columns $offset + 1 :
% offset + q(r)$ where $offset = \sum_{i=1}^{(r-1)}q(i) lev(i) + (k-1)*q(r)$.
% If idx = Gid{r} == k then
% Zs(idx,offset + 1 : offset + q(r)) = Z{r}(idx,:)
%
% Performing these steps on our dummy matrices produces the following Zs:
%
% Zs =
%
%      1     2     0     0    11    12    13     0     0     0     0     0     0
%      0     0     3     4     0     0     0    14    15    16     0     0     0
%      5     6     0     0     0     0     0     0     0     0    17    18    19
%      7     8     0     0     0     0     0     0     0     0    20    21    22
%      0     0     9    10     0     0     0    23    24    25     0     0     0
%
%
%
        
        % (1) Extract ZColNames, q, lev, GNames, Gid and GidLevelNames.
            ZColNames = model.RandomInfo.ZColNames;
                    q = model.RandomInfo.q;
                  lev = model.GroupingInfo.lev;
               GNames = model.GroupingInfo.GNames;
                  Gid = model.GroupingInfo.Gid;
        GidLevelNames = model.GroupingInfo.GidLevelNames;        
                 qlev = q.*lev;
                 
        % (2) Number of grouping variables.
        R = length(Gid);                         
                               
        % (3) Initialize a table to name the columns of Zs. The table
        % will have columns: Group, Level and Name. Here's a little example
        % of what this table looks like for grouping variable named 'Shape'
        % with 2 levels 'Circle' and 'Triangle'. Each level has 2
        % predictors named 'Intercept' and 'Slope'.
        %
        %   Group    Level     Name
        %   Shape    Circle    Intercept
        %   Shape    Circle    Slope
        %   Shape    Triangle  Intercept
        %   Shape    Triangle  Slope
        ZsColNames = table(cell(0,1),cell(0,1),cell(0,1),'VariableNames',{'Group','Level','Name'});
        
        % (4) Populate ZsColNames.
        for r = 1:R            
            groupname = repmat(GNames(r),qlev(r),1);
            levelname = repmat(GidLevelNames{r}',q(r),1);
            levelname = levelname(:);
            name = repmat(ZColNames{r}',lev(r),1);
            
            ZsColNames = [ZsColNames;table(groupname,levelname,name,...
                    'VariableNames',{'Group','Level','Name'})]; %#ok<AGROW>
        end
        
        % Debug: Ensure that ZsColNames using the new method are the same
        % as using the old method:
        %ZsColNamesOLD = makeSparseZNamesOLD(model);
        %assert( isequal(ZsColNames,ZsColNamesOLD) );
        
    end % end of makeSparseZNames.    

    function ZsColNames = makeSparseZNamesOLD(model)
%makeSparseZNames - Make names for the columns of big sparse Z matrix.        
%   ZsColNames = makeSparseZNames(model) takes an object model of type
%   LinearLikeMixedModel and makes ZsColNames, a table array such that each 
%   column of Zs returned by makeSparseZ is named by the corresponding row 
%   in ZsColNames. Every level of every grouping variable gets a random 
%   effects vector of the right size.
%
%   NOTE:
%
%   Z = R-by-1 cell array. Z{i} = N-by-q(i) design matrix for grouping variable i.
%   q = R-by-1 vector with q(i) = length of random effect vector for ith grouping variable.
% lev = R-by-1 array. lev(i) = number of levels in the ith grouping variable.
% Gid = R-by-1 cell array. Gid{i} = N-by-1 vector containing integers from 1 to lev(i)
%       indicating the membership of each observation into one of the levels from 1
%       to lev(i) of grouping variable i.
%
%%
% Suppose there are $r$ grouping variables. For the $i$ th grouping variable, we
% have $lev(i)$ levels. In each level of the $i$ th grouping variable, there is
% a random effect vector of length $q(i)$ which can be considered a draw from
% $N(0,\sigma^2 D_i)$. The co-variance $D_i$ is parameterized by a vector
% $\theta_i$ of length $lenTheta(i)$. The input vector $\theta$ is a concatenation
% of such vectors $\theta_i$. Thus:
%
% $\theta = [\theta_1; \theta_2; ...; \theta_r]$ and
%
% $length(\theta_i) = lenTheta(i)$
%
% If we concatenate all random effect vectors across all grouping variables
% into one big vector, we get the following vector:
%
% $b =
% [b(1)_1,b(1)_2,..,b(1)_{lev(1)},b(2)_1,b(2)_2,..,b(2)_{lev(2)}...,b(r)_1,b(r)_2,..,b(r)_{lev(r)}]^T$
%
% $b(i)_j = j$ th level random effect vector for $i$ th grouping variable
% of length $q(i) \times 1$
%
% The length of $b$ will be $\sum_{i = 1}^r q(i) lev(i)$.
%
% To understand, how Zs is constructed from Z, let us take a simple
% example. Suppose:
%
% Z{1} =
%
%     1     2
%     3     4
%     5     6
%     7     8
%     9    10
%
% Z{2} =
%
%    11    12    13
%    14    15    16
%    17    18    19
%    20    21    22
%    23    24    25
%
% lev(1) =
%
%     2
%
% lev(2) =
%
%     3
%
% Gid{1} =
%
%     1
%     2
%     1
%     1
%     2
%
% Gid{2} =
%
%     1
%     2
%     3
%     3
%     2
%
% The columns in matrix Zs correspond to the concatenated $b$ vector:
% $b =
% [b(1)_1,b(1)_2,..,b(1)_{lev(1)},b(2)_1,b(2)_2,..,b(2)_{lev(2)}...,b(r)_1,b(r)_2,..,b(r)_{lev(r)}]^T$
%
% All rows of Z{1} that satisfy Gid{1} == 1 go into columns 1:q(1) of Zs.
% If idx = Gid{1} == 1 then
% Zs(idx,1:q(1)) = Z{1}(idx,:)
%
% All rows of Z{r} that satisfy Gid{r} == k go into columns $offset + 1 :
% offset + q(r)$ where $offset = \sum_{i=1}^{(r-1)}q(i) lev(i) + (k-1)*q(r)$.
% If idx = Gid{r} == k then
% Zs(idx,offset + 1 : offset + q(r)) = Z{r}(idx,:)
%
% Performing these steps on our dummy matrices produces the following Zs:
%
% Zs =
%
%      1     2     0     0    11    12    13     0     0     0     0     0     0
%      0     0     3     4     0     0     0    14    15    16     0     0     0
%      5     6     0     0     0     0     0     0     0     0    17    18    19
%      7     8     0     0     0     0     0     0     0     0    20    21    22
%      0     0     9    10     0     0     0    23    24    25     0     0     0
%
%
%
        
        % (1) Extract ZColNames, q, lev, GNames, Gid and GidLevelNames.
            ZColNames = model.RandomInfo.ZColNames;
                    q = model.RandomInfo.q;
                  lev = model.GroupingInfo.lev;
               GNames = model.GroupingInfo.GNames;
                  Gid = model.GroupingInfo.Gid;
        GidLevelNames = model.GroupingInfo.GidLevelNames;        
        
        % (2) Number of grouping variables.
        R = length(Gid);                         
                               
        % (3) Initialize a table to name the columns of Zs. The table
        % will have columns: Group, Level and Name. Here's a little example
        % of what this table looks like for grouping variable named 'Shape'
        % with 2 levels 'Circle' and 'Triangle'. Each level has 2
        % predictors named 'Intercept' and 'Slope'.
        %
        %   Group    Level     Name
        %   Shape    Circle    Intercept
        %   Shape    Circle    Slope
        %   Shape    Triangle  Intercept
        %   Shape    Triangle  Slope
        ZsColNames = table(cell(0,1),cell(0,1),cell(0,1),'VariableNames',{'Group','Level','Name'});
        
        % (4) Populate Zs.
        for r = 1:R
            % For each grouping variable.
            for k = 1:lev(r)                                                
                % Name of grouping variable r in 1-by-1 cell GNames(r)
                groupname = cell(q(r),1);
                groupname(1:q(r)) = GNames(r);
                
                % Name of level k, grouping variable r: GidLevelNames{r}(k)
                levelname = cell(q(r),1);
                levelname(1:q(r)) = GidLevelNames{r}(k);
                
                % Predictors for Z{r}: ZColNames{r} of size 1-by-q(r)
                name = ZColNames{r}';
                
                % Append current naming info to ZsColnames
                ZsColNames = [ZsColNames;table(groupname,levelname,name,...
                    'VariableNames',{'Group','Level','Name'})]; %#ok<AGROW>
            end
        end                                                                                                                                                                                                                
        
    end % end of makeSparseZNamesOLD.    
    
    function X = fixedEffectsDesign(model)
%fixedEffectsDesign Returns the fixed effects design matrix.
%   X = fixedEffectsDesign(model) returns the N-by-P fixed effects design
%   matrix X for the linear like mixed effects model LME where N is the
%   number of observations and P is the number of fixed effects.

        % (1) Get subset of observations used in the fit and total number 
        % of observations N, whether used in the fit or not.
        subset = model.ObservationInfo.Subset;
        N = length(subset);
        
        % (2) How many fixed effects in the model.
        P = model.slme.p;

        % (3) Create a all NaN X of size N-by-P.
        X = NaN(N,P);
        
        % (4) Now fill in info for observations used in the fit.        
        X(subset,:) = model.FixedInfo.X;

    end % end of fixedEffectsDesign.    

    function [Z,gnames] = randomEffectsDesign(model,gnumbers)
%randomEffectsDesign Extract random effects design matrices.
%   Z = randomEffectsDesign(LME) returns the overall random effects design
%   matrix corresponding to a vector B of all random effects in the linear
%   mixed effects model LME. Suppose LME has R grouping variables named
%   g_1,...,g_R. Let Q_1,...,Q_R be the length of random effects vectors
%   associated with g_1,...,g_R respectively. Also, suppose g_1,...,g_R
%   have levels M_1,...,M_R respectively. Then B will be a column vector of
%   length Q_1*M_1 + Q_2*M_2 + ... + Q_R*M_R. B is made by concatenating
%   the BLUPs of random effects vectors corresponding to each level of each
%   grouping variable in the following order:
%
%   g_1 level 1, g_1 level 2,..., g_1 level M_1 followed by
%   g_2 level 1, g_2 level 2,..., g_2 level M_2 followed by
%      .            .                .           
%   g_R level 1, g_R level 2,..., g_R level M_R
%
%   If LME has N observations then Z will be of size N-by-length(B) such
%   that Z * B is a N-by-1 vector that represents the contribution of all
%   random effects to the response of the LME.
%
%   ZSUB = designMatrix(LME,GNUMBERS) returns a submatrix of the full
%   random effects design matrix. GNUMBERS is a length K integer array with
%   elements between 1 and R. ZSUB is a subset of the full random effects
%   design matrix corresponding to the grouping variable names indicated by
%   integers in GNUMBERS. For example, suppose GNUMBERS is [1,R] then this
%   specifies only grouping variables g_1 and g_R. Let BSUB be a vector
%   made by concatenating BLUPs of random effects vectors corresponding to
%   each level of g_1 and g_R in the following order:
%
%   g_1 level 1, g_1 level 2,..., g_1 level M_1 followed by         
%   g_R level 1, g_R level 2,..., g_R level M_R
%
%   Then DSUB will be a N-by-length(BSUB) matrix such that DSUB*BSUB
%   represents the contribution of all random effects corresponding to
%   grouping variables g_1 and g_R to the response of the LME. If GNUMBERS
%   is empty, the full random effects design matrix is returned.
%
%   [ZSUB,GNAMES] = designMatrix(LME,GNUMBERS) also returns a K-by-1 cell
%   array containing the names of grouping variables corresponding to
%   integers in GNUMBERS.
%
%   This method also applies to a GeneralizedLinearMixedModel object.

        % (1) If gnumbers is not given, set gnumbers to []. This marks that
        % we want the full random effects design matrix.
        if nargin == 1
            gnumbers = [];
        end

        % (2) How many groups do we have.
        R = model.GroupingInfo.R;
        
        % (3) Get subset of observations used in the fit and total number 
        % of observations N, whether used in the fit or not.
        subset = model.ObservationInfo.Subset;
        N = length(subset);               
                
        if isempty(gnumbers)
            % (4) If gnumbers is [], set gnumbers equal to 1:R where R is
            % the total number of grouping variables.
            gnumbers = (1:R)';
        else
            % (5) Validate gnumbers.
            gnumbers = model.validateGNumbers(gnumbers,R);        
        end
                        
        % (6) If R = 0, we are done.
        if R == 0
            Z = sparse(N,0);
            gnames = {};
            return;
        end        
        
        % (7) Extract lev and q.
        lev = model.GroupingInfo.lev;
          q = model.RandomInfo.q;

        % (8) Form the product q * lev (element by element).
        qlev = q .* lev;
                      
        % (9) Loop over integers in gnumbers and build an index vector idx 
        % which will be used to extract relevant columns from model.slme.Z.       
        idx = [];
        gnames = cell(length(gnumbers), 1);
        for r = 1:length(gnumbers)
            
            % (1) Get index of grouping variable.
            k = gnumbers(r);
            
            % (2) gnames{r} is the name associated with grouping variable
            % with index k.
            gnames{r} = model.GroupingInfo.GNames{k};

            % (3) Grouping variable with index k is contained in the
            % following columns of model.slme.Z:
            %
            %   cols = offset + 1 : offset + q(k)*lev(k) where
            % offset = sum_{i=1}^(k-1) q(i)*lev(i)
            
            offset = sum(qlev(1:(k-1)));
            cols = offset + 1 : offset + qlev(k);
            
            % (4) Append the cols vector to the overall index vector idx
            % for all grouping variables.
            idx = [idx,cols];             %#ok<AGROW>      
            
        end
        
        % (10) Now extract columns from model.slme.Z corresponding to the 
        % index vector idx.
        Zr = model.slme.Z(:,idx);                
               
        % (11) Create a all zero N-by-size(Zr,2) sparse matrix Z.
        Z = sparse(N,size(Zr,2));
        
        % (12) Now fill in info for observations used in the fit.        
        Z(subset,:) = Zr;        
        
    end % end of randomEffectsDesign.    

    function [X,XColNames,XCols2Terms] = createDesignMatrix(model,ds,F)
%createDesignMatrix - Create a design matrix from a dataset/table array.
%   [X,XColNames,XCols2Terms] = createDesignMatrix(model,ds,F) takes a
%   dataset/table ds, a LinearFormula F and returns the N-by-P design
%   matrix X, a cell array XColNames of length P containing the name of
%   each column in X and an integer array XCols2Terms of size 1-by-P
%   indicating which "term" a particular column in X belongs to. All
%   columns of X that belong to the same "term" can be entered into an
%   ANOVA table with an F-test.
        
        % (1) Ensure that ds has all predictor variables. ds may not
        % contain all variables in table model.Variables. But ds must
        % contain all variable names in model.PredictorNames.
        if isa(ds,'dataset')
            ds = dataset2table(ds);
        end
          
        tf = ismember(model.PredictorNames,ds.Properties.VariableNames);        
        if ~all(tf)
            error(message('stats:classreg:regr:TermsRegression:MissingVariable'));
        end
        
        % (2) Get column indices of ds VarNames in model VariableNames.        
        [tf,varLocs] = ismember(ds.Properties.VariableNames,model.VariableNames);
        if ~all(tf)
            % Remove extra variables from ds so that we don't have to
            % construct Terms, IsCategorical, Range, and DummyVarCoding
            % values.
            ds = ds(:,tf);
        end
        varLocs = varLocs(tf);

        % (3) Get terms matrix specified by F.
        terms = F.Terms;
        
        % (4) Ensure that F.VariableNames matches ds.VariableNames.
        assert( isequal(F.VariableNames,model.VariableNames) | ...
            isequal(F.VariableNames,model.VariableNames') );
        
        % (5) Create design matrix using terms and varLocs.
        [X,~,~,XCols2Terms,XColNames] ...
            = classreg.regr.modelutils.designmatrix(ds,'Model',terms(:,varLocs), ...
            'DummyVarCoding',model.DummyVarCoding, ...
            'CategoricalVars',model.VariableInfo.IsCategorical(varLocs), ...
            'CategoricalLevels',model.VariableInfo.Range(varLocs));                        
        
    end % end of createDesignMatrix.
    
end

% Abstract utility methods.
methods (Abstract,Access={?classreg.regr.LinearLikeMixedModel})
   
    np = getTotalNumberOfParameters(model);
    
end

% Shared static utility methods.
methods(Static,Access='protected')
    
    function [G,GName] = makeInteractionVar(ds,interactionVars)
%makeInteractionVar - Makes an interaction grouping variable.
%   [G,GName] = makeInteractionVar(ds,interactionVars) takes a dataset/table 
%   array ds containing all variables listed in a cell array of strings 
%   interactionVars and returns a single nominal variable G that represents
%   the interaction between variables in interactionVars. GName is the name 
%   assigned to G. GName is created by putting a ':' between the strings 
%   in interactionVars. If interactionVars is {} then both G and GName are
%   [].

        if isa(ds,'dataset')
            ds = dataset2table(ds);
        end
        
        % (1) Ensure that interactionVars is a cell array of strings.
        assert( iscellstr(interactionVars) );
        
        % (2) Ensure that all variables in interactionVars also appear in
        % the table ds.
        assert( all(ismember(interactionVars,ds.Properties.VariableNames)) );
       
        % (3) Get the required interaction by first converting to nominal
        % and then taking the .* product on nominal vectors.
        k = length(interactionVars);
        if k >= 1
            G = nominal(ds.(interactionVars{1}));
            GName = interactionVars{1};
            for i = 2:k
                G = G.*nominal(ds.(interactionVars{i}));
                GName = [GName,':',interactionVars{i}];         %#ok<AGROW>
            end
        else
            G = [];
            GName = [];
        end
        
        % (4) Drop unused levels from G.
        G = removecats(G);
        
    end % end of makeInteractionVar.    

    function newgid = reorderGroupIDs(gid,names,newnames)
%reorderGroupIDs - Reorders a group id vector.       
%   newgid = reorderGroupIDs(model,gid,names,newnames) takes a N-by-1
%   integer vector gid and a M-by-1 cell array names.
%
%   (1) Levels with gid == j have name names{j}.
%
%   (2) If names{j} is equal to newnames{k} then levels with gid == j are
%   assigned the id k in newgid.
%
%   (3) If names{j} is not found in newnames{k} then levels with gid == j
%   are assigned the id NaN in newgid.
%
%   newgid has the same length as gid.

        % (1) gid must be an integer vector.
        assert( internal.stats.isIntegerVals(gid,1,Inf) );
        assert( isvector(gid) );
        if size(gid,1) == 1
            gid = gid';
        end
        
        % (2) names and newnames must be cell vector of strings.
        assert( iscellstr(names) & isvector(names) );
        if size(names,1) == 1
            names = names';
        end
        assert( iscellstr(newnames) & isvector(newnames) );
        if size(newnames,1) == 1
            newnames = newnames';
        end
        
        % (3) Consider [tf,newintid] = ismember( names, newnames );
        % tf(j) will be true if names{j} is found in newnames. newintid(j)
        % will contain the location k where names{j} was found in newnames
        % or 0 if names{j} was not found in newnames.
        [~,newintid] = ismember( names, newnames );

        % (4) Form the output vector newgid.
        newgid = NaN(size(gid));        
        for j = 1:length(newintid)
            % index k such that names{j} matches newnames{k}.
            k = newintid(j);
            if k ~= 0      
                newgid( gid == j ) = k;                
            end           
        end
 
    end % end of reorderGroupIDs.
        
    function tf = isMatrixNested(Xsmall,Xbig)
%isMatrixNested - Checks if the span of Xbig contains Xsmall.        
%   tf = isMatrixNested(Xsmall,Xbig) takes two matrices Xsmall and Xbig 
%   with the same number of rows and checks if the span of Xbig contains
%   Xsmall. If so, tf is true and if not, tf is false.

        % (1) Make sure that Xsmall and Xbig are numeric, real matrices.
        assert( isnumeric(Xsmall) & isreal(Xsmall) & ismatrix(Xsmall) );
        assert( isnumeric(Xbig)   & isreal(Xbig)   & ismatrix(Xbig)   );
        
        % (2) Ensure that Xsmall and Xbig have the same number of rows.
        assert( size(Xsmall,1) == size(Xbig,1) );

        % (3) Suppose the span of Xbig contains Xsmall. Then we can write
        % Xbig * G = Xsmall for some matrix G. We would like to find an
        % estimate GHat for G in the "least squares" sense and then form an
        % estimate of Xsmall using Xsmall_recon = Xbig * GHat.
        
            % 3.1 If Xbig is a square matrix, add a row of zeros to Xbig
            % and Xsmall to force least squares solution in mldivide.
            [Nbig,qbig] = size(Xbig);
            if ( Nbig == qbig )
                Xbig(Nbig+1,:)   = 0;
                Xsmall(Nbig+1,:) = 0;
            end
        
            % 3.2 Form the estimate GHat.
            warnState = warning('query','all');
            warning('off','MATLAB:nearlySingularMatrix');
            warning('off','MATLAB:illConditionedMatrix');
            warning('off','MATLAB:singularMatrix');
            warning('off','MATLAB:rankDeficientMatrix');
            cleanupObj = onCleanup(@() warning(warnState));
            GHat = Xbig \ Xsmall;
            
            % 3.3 Form the reconstruction Xsmall_recon.
            Xsmall_recon = Xbig*GHat;
            
        % (4) Measure the error between Xsmall and Xsmall_recon.
        err = max(max(abs(Xsmall - Xsmall_recon)));

        % (5) err should be close to 0 or [] (happens when Xsmall is []).
        if ( isempty(err) || err <= sqrt(eps(class(Xbig)))*max(size(Xbig)) )
            tf = true;
        else
            tf = false;
        end                    
        
    end % end of isMatrixNested.
    
    function tf = isMatrixNestedOLD(Xsmall,Xbig)
%isMatrixNested - Checks if the span of Xbig contains Xsmall.        
%   tf = isMatrixNested(Xsmall,Xbig) takes two matrices Xsmall and Xbig 
%   with the same number of rows and checks if the span of Xbig contains
%   Xsmall. If so, tf is true and if not, tf is false.

        % (1) Make sure that Xsmall and Xbig are numeric, real matrices.
        assert( isnumeric(Xsmall) & isreal(Xsmall) & ismatrix(Xsmall) );
        assert( isnumeric(Xbig)   & isreal(Xbig)   & ismatrix(Xbig)   );
        
        % (2) Ensure that Xsmall and Xbig have the same number of rows.
        assert( size(Xsmall,1) == size(Xbig,1) );

        % (3) Convert Xbig to full matrix (so that orth can work).
        Xbig = full(Xbig);
        
        % (4) Find an orthogonal basis for the columns of Xbig such that
        % Qbig'*Qbig = I (identity).
        Qbig = orth(Xbig);
        
        % (5) If Xsmall = Qbig * G for some G holds then G = Qbig'*Xsmall.
        % So Xsmall = Qbig * (Qbig'*Xsmall) must hold.
        Xsmall_recon = Qbig*(Qbig'*Xsmall);
        
        % (6) Ensure Xsmall_recon is numerically equal to Xsmall.
        err = max(abs( Xsmall(:) - Xsmall_recon(:) ));
        
        % (7) err should be close to 0 or [] (happens when Xsmall is []).
        if ( isempty(err) || err <= sqrt(eps(class(Xbig)))*max(size(Xbig)) )
            tf = true;
        else
            tf = false;
        end        
        
    end % end of isMatrixNestedOLD.

    function lrt = standardLRT(smallModel,bigModel,smallModelName,bigModelName)
%standardLRT - Standard likelihood ratio test.        
%   lrt = standardLRT(smallModel,bigModel,smallModelName,bigModelName)      
%   takes LinearLikeMixedModel objects smallModel and bigModel with names
%   smallModelName and bigModelName respectively. smallModel is assumed to
%   be nested in bigModel. The output lrtable is a dataset with the
%   following columns:
%
%           Model       The name of the model.
%           DF          The number of free parameters in the model.
%           AIC         Akaike information criterion for the model.
%           BIC         Bayesian information criterion for the model.
%           LogLik      The maximized log-likelihood for the model.
%           LRStat      Likelihood ratio test statistic for comparing 
%                       bigModel versus smallModel.
%           deltaDF     DF for bigModel minus DF for smallModel.
%           pValue      p-value for the likelihood ratio test.        

% NOTE: The ModelCriterion property looks like this:
%
%           AIC       BIC       LogLikelihood    Deviance
%          841.21    861.97    -412.61          825.21  
%
% It may be a structure, instead of a table but it will always have
% properties AIC, BIC, LogLikelihood and Deviance.
%
%
% The output lrt should be a table that looks like this:
%
% Model     DF   AIC     BIC     LogLik    LRStat  deltaDF   pValue
% lme       4    149.22  156.17  -70.609
% altlme    6    149.43  159.85  -68.714   3.7896    2       0.15035   
%

        % (1) Initialize the output table.
        lrt = table();
        
        % (2) Add the Model column.
        modelName = cell(2,1);
        modelName{1} = smallModelName;
        modelName{2} = bigModelName;
        lrt.Model = nominal(modelName);

        % (3) Add the DF column.
        DF = zeros(2,1);
        DF(1) = getTotalNumberOfParameters(smallModel);
        DF(2) = getTotalNumberOfParameters(bigModel);
        lrt.DF = DF;
        
        % (3) Get the model criterion values for all models.
        smallModelCrit = smallModel.ModelCriterion;
          bigModelCrit =   bigModel.ModelCriterion;
        
        % (4) Add the AIC column.
        AIC = zeros(2,1);
        AIC(1) = smallModelCrit.AIC;
        AIC(2) =   bigModelCrit.AIC;
        lrt.AIC = AIC;
        
        % (5) Add the BIC column.
        BIC = zeros(2,1);
        BIC(1) = smallModelCrit.BIC;
        BIC(2) =   bigModelCrit.BIC;
        lrt.BIC = BIC;
        
        % (6) Add the LogLik column.
        LogLik = zeros(2,1);
        LogLik(1) = smallModelCrit.LogLikelihood;
        LogLik(2) =   bigModelCrit.LogLikelihood;
        lrt.LogLik = LogLik;
        
        % (7) Add the LRStat column.
        LRatio = zeros(2,1);
        LRatio(1) = 0;
        LRatio(2) = 2*(LogLik(2) - LogLik(1));
        LRatioAbsent = [true;false];
        lrt.LRStat = internal.stats.DoubleTableColumn(LRatio,LRatioAbsent);
        
        % (8) Add the deltaDF column.
        deltaDF = zeros(2,1);
        deltaDF(1) = 0;
        deltaDF(2) = DF(2) - DF(1);
        deltaDFAbsent = [true;false];
        lrt.deltaDF = internal.stats.DoubleTableColumn(deltaDF,deltaDFAbsent);
        
        % (9) Add the pValue column.
        pValue = zeros(2,1);
        pValue(1) = 0;
        pValue(2) = 1 - chi2cdf( LRatio(2), deltaDF(2) );
        pValueAbsent = [true;false];
        lrt.pValue = internal.stats.DoubleTableColumn(pValue,pValueAbsent);

        % (10) Add a title to lrt.
        ttl = getString(message('stats:LinearMixedModel:Title_LRT'));
        lrt = classreg.regr.lmeutils.titleddataset(lrt,ttl);
               
    end % end of standardLRT.
    
    function tds = removeTitle(tds)
%removeTitle - Removes the title from a titled dataset.        
%   tds = removeTitle(tds) takes a titled dataset tds and removes title. If
%   tds is a titled dataset then the output is also a titled dataset with
%   the title removed. If tds is not a titled dataset then the output is
%   identical to the input.

        if isa(tds,'classreg.regr.lmeutils.titleddataset')
            tds = settitle(tds,'');
        end

    end % end of convertTitledDatasetToDataset.
    
    function str = formatBold(str)
%formatBold - Format the string str so that it gets displayed in bold.
%   str = formatBold(str) adds a <strong> and </strong> tag around the 
%   string str if feature('hotlinks') is true.       

        if feature('hotlinks')
            str = ['<strong>',str,'</strong>'];
        end
        
    end % end of formatBold.    
    
    function [haveDataset,ds,X,Z,G,otherArgs] = handleDatasetOrMatrixInput(varargin)
%handleDatasetOrMatrixInput - Handle dataset/table or matrix input.       
%   [haveDataset,ds,X,Z,G,otherArgs] = handleDatasetOrMatrixInput(varargin)
%   takes a variable number of input arguments and returns info on the
%   dataset/table or matrix input. The random and predict methods may be called
%   as random(lme,DS) or random(lme,X) or random(lme,X,Z) or
%   random(lme,X,Z,G). This is a convenience method to figure out which
%   case we are dealing with. Here's what's returned:
%
%   haveDataset - true if input was a dataset/table.
%            ds - input dataset/table if input was a dataset/table or [].
%             X - matrix input X if supplied or [] if not supplied.
%             Z - matrix input Z if supplied or [] if not supplied.
%             G - matrix input G if supplied or [] if not supplied.
%     otherArgs - other input args.
        

        % (1) length of varargin must be >= 1.
        assert( length(varargin) >= 1 );

        % (2) Do we have a dataset/table in the input?
        if isa(varargin{1},'dataset')
            varargin{1}=dataset2table(varargin{1});
        end
        haveDataset = isa(varargin{1},'table');
        
        % (3) Handle the various cases.
        if haveDataset == true
            ds = varargin{1};
            X = [];
            Z = [];
            G = [];
            otherArgs = varargin(2:end);
        else
            nargs = length(varargin);            
            switch nargs
                case 1
                    % Get X, set ds, Z, G = [].
                    ds = [];
                    X = varargin{1};
                    Z = [];
                    G = [];
                    otherArgs = varargin(2:end);
                case 2
                    % Get X, Z, set ds, G = [].
                    ds = [];
                    X = varargin{1};
                    Z = varargin{2};
                    G = [];
                    otherArgs = varargin(3:end);
                case 3
                    % Get X, Z and G, set ds = [].
                    ds = [];
                    X = varargin{1};
                    Z = varargin{2};
                    G = varargin{3};
                    otherArgs = varargin(4:end);
                otherwise
                    % Get X and Z, set ds = [].
                    ds = [];
                    X = varargin{1};
                    Z = varargin{2};
                    % Try to get G. If varargin{3} is a string, and the 
                    % length of rest of the args is odd then it must be 
                    % a parameter name/value item, otherwise varargin{3} is
                    % a potential G.                    
                    if internal.stats.isString(varargin{3}) && rem(length(varargin(4:end)),2) == 1
                        G = [];                 
                        otherArgs = varargin(3:end);
                    else                        
                        G = varargin{3};
                        otherArgs = varargin(4:end);
                    end                       
            end                
        end        
        
    end % end of handleDatasetOrMatrixInput.
    
end

% Shared static utility methods (public, hidden)
methods(Static,Access='public',Hidden=true)
    
    function Zs = makeSparseZ(Z,q,lev,Gid,N)
%makeSparseZ - Make big sparse Z matrix for LME/GLME in standard form.        
%   Zs = makeSparseZ(Z,q,lev,Gid,N) takes an object model of type
%   LinearMixedModel and makes a big sparse matrix Zs that is used for
%   fitting the LME in standard form. Every level of every grouping
%   variable gets a random effects vector of the right size.
%
%
%   NOTE:
%
%   Z = R-by-1 cell array. Z{i} = N-by-q(i) design matrix for grouping variable i.
%   q = R-by-1 vector with q(i) = length of random effect vector for ith grouping variable.
% lev = R-by-1 array. lev(i) = number of levels in the ith grouping variable.
% Gid = R-by-1 cell array. Gid{i} = N-by-1 vector containing integers from 1 to lev(i)
%       indicating the membership of each observation into one of the levels from 1
%       to lev(i) of grouping variable i. Gid{i} is such that integer j in Gid{i}
%       has name model.GroupingInfo.GidLevelNames{i}{j}.
%   N = number of observations to be used in the fit.
%%
% Suppose there are $r$ grouping variables. For the $i$ th grouping variable, we
% have $lev(i)$ levels. In each level of the $i$ th grouping variable, there is
% a random effect vector of length $q(i)$ which can be considered a draw from
% $N(0,\sigma^2 D_i)$. The co-variance $D_i$ is parameterized by a vector
% $\theta_i$ of length $lenTheta(i)$. The input vector $\theta$ is a concatenation
% of such vectors $\theta_i$. Thus:
%
% $\theta = [\theta_1; \theta_2; ...; \theta_r]$ and
%
% $length(\theta_i) = lenTheta(i)$
%
% If we concatenate all random effect vectors across all grouping variables
% into one big vector, we get the following vector:
%
% $b =
% [b(1)_1,b(1)_2,..,b(1)_{lev(1)},b(2)_1,b(2)_2,..,b(2)_{lev(2)}...,b(r)_1,b(r)_2,..,b(r)_{lev(r)}]^T$
%
% $b(i)_j = j$ th level random effect vector for $i$ th grouping variable
% of length $q(i) \times 1$
%
% The length of $b$ will be $\sum_{i = 1}^r q(i) lev(i)$.
%
% To understand, how Zs is constructed from Z, let us take a simple
% example. Suppose:
%
% Z{1} =
%
%     1     2
%     3     4
%     5     6
%     7     8
%     9    10
%
% Z{2} =
%
%    11    12    13
%    14    15    16
%    17    18    19
%    20    21    22
%    23    24    25
%
% lev(1) =
%
%     2
%
% lev(2) =
%
%     3
%
% Gid{1} =
%
%     1
%     2
%     1
%     1
%     2
%
% Gid{2} =
%
%     1
%     2
%     3
%     3
%     2
%
% The columns in matrix Zs correspond to the concatenated $b$ vector:
% $b =
% [b(1)_1,b(1)_2,..,b(1)_{lev(1)},b(2)_1,b(2)_2,..,b(2)_{lev(2)}...,b(r)_1,b(r)_2,..,b(r)_{lev(r)}]^T$
%
% All rows of Z{1} that satisfy Gid{1} == 1 go into columns 1:q(1) of Zs.
% If idx = Gid{1} == 1 then
% Zs(idx,1:q(1)) = Z{1}(idx,:)
%
% All rows of Z{r} that satisfy Gid{r} == k go into columns $offset + 1 :
% offset + q(r)$ where $offset = \sum_{i=1}^{(r-1)}q(i) lev(i) + (k-1)*q(r)$.
% If idx = Gid{r} == k then
% Zs(idx,offset + 1 : offset + q(r)) = Z{r}(idx,:)
%
% Performing these steps on our dummy matrices produces the following Zs:
%
% Zs =
%
%      1     2     0     0    11    12    13     0     0     0     0     0     0
%      0     0     3     4     0     0     0    14    15    16     0     0     0
%      5     6     0     0     0     0     0     0     0     0    17    18    19
%      7     8     0     0     0     0     0     0     0     0    20    21    22
%      0     0     9    10     0     0     0    23    24    25     0     0     0
%
%        
%                     Z = model.RandomInfo.Z;
%                     q = model.RandomInfo.q;
%                   lev = model.GroupingInfo.lev;
%                   Gid = model.GroupingInfo.Gid;

        % 1. Form the product q * lev (element by element).
        qlev = q .* lev;
        
        % 2. Number of grouping variables. Return early if no grouping
        % variables.
        R = length(Gid);
        if (R == 0)
            Zs = sparse(N,0);
            return;
        end
                        
        % 3. Populate Zs. Zs is a sparse matrix of size N-by-sum(qlev).
        % Accumulate rowidx, colidx, validx for later call to sparse.
        rowidx = cell(sum(lev),1);
        colidx = rowidx;
        validx = rowidx;
        
            count = 1;
            for r = 1:R
                % For each grouping variable.        

                % 3.1 If Gid{r} contains any NaN's, replace them by
                % lev(r)+1. We do this so that subsequent call to
                % accumarray works.
                Gid{r}(isnan(Gid{r})) = lev(r) + 1;

                % 3.2 Create a cell array idxr such that:
                % 
                % idxr{k} = find(Gid{r}==k) 
                %
                % using accumarray. We could call accumarray like this:
                %
                % accumarray(Gid{r},(1:N)',[],@(x) {})
                %
                % Alternatively, we can sort Gid{r} first and apply the
                % sorting to (1:N)' before calling accumarray.
                vals     = (1:N)';
                [subs,I] = sortrows(Gid{r});
                idxr     = accumarray(subs,vals(I),[],@(x) {x});
                lenidxr  = length(idxr);
                if lenidxr < lev(r)
                    % Pad so that length(idxr) is equal to lev(r).
                    idxr(lenidxr+1:lev(r)) = {zeros(0,1)};
                end

                for k = 1:lev(r)                
                    % For each level of grouping variable r.

                    idx = idxr{k};

                    if isempty(idx)
                        rowidx{count} = zeros(0,1);
                        colidx{count} = zeros(0,1);
                        validx{count} = zeros(0,1);
                    else                
                        offset = sum(qlev(1:(r-1))) + (k-1)*q(r);                

                        % 3.3 Here's what we want to do:
                        %
                        % Zs(idx, offset + 1 : offset + q(r)) = Z{r}(idx,:);
                        %
                        % Since sparse indexing is slow, our idea is to use
                        % the syntax S = sparse(i,j,s,m,n,nzmax) to create
                        % a sparse matrix where i, j and s are vectors such
                        % that S is created by: S(i(k),j(k)) = s(k).
                        %
                        % Here's an example:
                        %
                        % A = [1     2     3     4     5
                        %      6     7     8     9    10
                        %     11    12    13    14    15];
                        %
                        % row = [1
                        %        2];
                        %
                        % col = [2,3,4];
                        %
                        % val = [101   102   103
                        %        104   105   106]
                        %
                        % Clearly, A(row,col) = val produces:
                        %
                        % A = [1   101   102   103     5
                        %      6   104   105   106    10
                        %     11    12    13    14    15];
                        % 
                        % Consider the commands:
                        %
                        % nrow = length(row);
                        % ncol = length(col);
                        % row  = repmat(row,1,ncol);
                        % col  = repmat(col,nrow,1);
                        % 
                        % then row and col look like this:
                        %
                        % row = [1     1     1
                        %        2     2     2];
                        %
                        % col = [2     3     4
                        %        2     3     4];
                        %
                        % Finally, [row(:),col(:),val(:)] gives:
                        %
                        % [row(:),col(:),val(:)] = 
                        %
                        % [1     2   101
                        %  2     2   104
                        %  1     3   102
                        %  2     3   105
                        %  1     4   103
                        %  2     4   106];
                        %
                        % Each row of the above output correctly produces
                        % the row index, column index and value
                        % combination. We can loop over r = 1:R, k =
                        % 1:lev(r) and collect vectors row(:) into rowidx,
                        % col(:) into colidx and val(:) into validx.
                        % Finally, we can create Zs like this:
                        %
                        % Zs = sparse(rowidx,colidx,validx,N,sum(qlev));                

                        row = idx;
                        col = offset + 1 : offset + q(r);
                        nrow = length(row);
                        ncol = length(col);

                        row = repmat(row,1,ncol);
                        col = repmat(col,nrow,1);
                        val = Z{r}(idx,:);

                        row = row(:);
                        col = col(:);
                        val = val(:);

                        rowidx{count} = row;
                        colidx{count} = col;
                        validx{count} = val;                               
                    end

                    count = count + 1;

                end
            end                                                                                                                                                                                                                        
        
        % 4. Make the big sparse Zs matrix.
        rowidx = cell2mat(rowidx);
        colidx = cell2mat(colidx);
        validx = cell2mat(validx);        
        Zs = sparse(rowidx,colidx,validx,N,sum(qlev));
        
        % Debug: Ensure that the new way of making Zs gives the same
        % results as the old way.
        %ZsOLD = classreg.regr.LinearLikeMixedModel.makeSparseZOLD(Z,q,lev,Gid,N);
        %assert( isequal(Zs,ZsOLD) );
        
    end % end of makeSparseZ.

    function Zs = makeSparseZOLD(Z,q,lev,Gid,N)
%makeSparseZ - Make big sparse Z matrix for LME/GLME in standard form.        
%   Zs = makeSparseZ(Z,q,lev,Gid,N) takes an object model of type
%   LinearMixedModel and makes a big sparse matrix Zs that is used for
%   fitting the LME in standard form. Every level of every grouping
%   variable gets a random effects vector of the right size.
%
%
%   NOTE:
%
%   Z = R-by-1 cell array. Z{i} = N-by-q(i) design matrix for grouping variable i.
%   q = R-by-1 vector with q(i) = length of random effect vector for ith grouping variable.
% lev = R-by-1 array. lev(i) = number of levels in the ith grouping variable.
% Gid = R-by-1 cell array. Gid{i} = N-by-1 vector containing integers from 1 to lev(i)
%       indicating the membership of each observation into one of the levels from 1
%       to lev(i) of grouping variable i. Gid{i} is such that integer j in Gid{i}
%       has name model.GroupingInfo.GidLevelNames{i}{j}.
%   N = number of observations to be used in the fit.
%%
% Suppose there are $r$ grouping variables. For the $i$ th grouping variable, we
% have $lev(i)$ levels. In each level of the $i$ th grouping variable, there is
% a random effect vector of length $q(i)$ which can be considered a draw from
% $N(0,\sigma^2 D_i)$. The co-variance $D_i$ is parameterized by a vector
% $\theta_i$ of length $lenTheta(i)$. The input vector $\theta$ is a concatenation
% of such vectors $\theta_i$. Thus:
%
% $\theta = [\theta_1; \theta_2; ...; \theta_r]$ and
%
% $length(\theta_i) = lenTheta(i)$
%
% If we concatenate all random effect vectors across all grouping variables
% into one big vector, we get the following vector:
%
% $b =
% [b(1)_1,b(1)_2,..,b(1)_{lev(1)},b(2)_1,b(2)_2,..,b(2)_{lev(2)}...,b(r)_1,b(r)_2,..,b(r)_{lev(r)}]^T$
%
% $b(i)_j = j$ th level random effect vector for $i$ th grouping variable
% of length $q(i) \times 1$
%
% The length of $b$ will be $\sum_{i = 1}^r q(i) lev(i)$.
%
% To understand, how Zs is constructed from Z, let us take a simple
% example. Suppose:
%
% Z{1} =
%
%     1     2
%     3     4
%     5     6
%     7     8
%     9    10
%
% Z{2} =
%
%    11    12    13
%    14    15    16
%    17    18    19
%    20    21    22
%    23    24    25
%
% lev(1) =
%
%     2
%
% lev(2) =
%
%     3
%
% Gid{1} =
%
%     1
%     2
%     1
%     1
%     2
%
% Gid{2} =
%
%     1
%     2
%     3
%     3
%     2
%
% The columns in matrix Zs correspond to the concatenated $b$ vector:
% $b =
% [b(1)_1,b(1)_2,..,b(1)_{lev(1)},b(2)_1,b(2)_2,..,b(2)_{lev(2)}...,b(r)_1,b(r)_2,..,b(r)_{lev(r)}]^T$
%
% All rows of Z{1} that satisfy Gid{1} == 1 go into columns 1:q(1) of Zs.
% If idx = Gid{1} == 1 then
% Zs(idx,1:q(1)) = Z{1}(idx,:)
%
% All rows of Z{r} that satisfy Gid{r} == k go into columns $offset + 1 :
% offset + q(r)$ where $offset = \sum_{i=1}^{(r-1)}q(i) lev(i) + (k-1)*q(r)$.
% If idx = Gid{r} == k then
% Zs(idx,offset + 1 : offset + q(r)) = Z{r}(idx,:)
%
% Performing these steps on our dummy matrices produces the following Zs:
%
% Zs =
%
%      1     2     0     0    11    12    13     0     0     0     0     0     0
%      0     0     3     4     0     0     0    14    15    16     0     0     0
%      5     6     0     0     0     0     0     0     0     0    17    18    19
%      7     8     0     0     0     0     0     0     0     0    20    21    22
%      0     0     9    10     0     0     0    23    24    25     0     0     0
%
%
%
        
        % (1) Extract Z, q, lev, and Gid.
%                     Z = model.RandomInfo.Z;
%                     q = model.RandomInfo.q;
%                   lev = model.GroupingInfo.lev;
%                   Gid = model.GroupingInfo.Gid;

        % (2) Form the product q * lev (element by element).
        qlev = q .* lev;
        
        % (3) Number of grouping variables.
        R = length(Gid);
        
        % (4) Number of observations.
        %N = model.NumObservations;
        %N = size(Z{1},1);
        
        % (5) Initialize Zs - a sparse matrix.
        Zs = sparse(N,sum(qlev));
                        
        % (6) Populate Zs.
        for r = 1:R
            % For each grouping variable.
            for k = 1:lev(r)
                % For each level of grouping variable r.
                idx = ( Gid{r} == k );
                offset = sum(qlev(1:(r-1))) + (k-1)*q(r);
                Zs(idx, offset + 1 : offset + q(r)) = Z{r}(idx,:); %#ok<SPRIX>                                                
            end
        end                                                                                                                                                                                                                        
        
    end % end of makeSparseZOLD.
    
end

% Shared static utility methods for validation.
methods(Static,Access='protected')  
    
    function covariancepattern = ...
            validateCovariancePattern(covariancepattern,R,Q)
%validateCovariancePattern - Validates the 'CovariancePattern' option.
%   covariancepattern = validateCovariancePattern(covariancepattern,R,Q)
%   takes a tentative value for 'CovariancePattern' option, the number of
%   groups R, and Q, a length R vector containing the size of random 
%   effects vectors in each group and either returns the validated 
%   covariancepattern (a cell array of length R-by-1) or throws an error 
%   message.
%
%   What is checked?
%   
%   Either:
%
%   (1) If R = 1 then covariancepattern is a string contained in the list 
%   LinearMixedModel.AllowedCovariancePatterns or a Q(1)-by-Q(1) logical 
%   matrix.
%
%       or
%
%   (2) If R >= 1 then covariancepattern is a length R cell array such that 
%   covariancepattern{i} is either a string contained in the list 
%   LinearMixedModel.AllowedCovariancePatterns or a Q(i)-by-Q(i) logical
%   matrix.
%
%   If both (1) and (2) are not satisfied, then an error message is thrown.

        % (1) Ensure that R >= 0 and R is an integer.
        assert(R >= 0 & internal.stats.isScalarInt(R));
        
        % (2) Ensure that Q is an integer vector with elements in the range
        % 0 to Inf. Also length(Q) must be equal to R.
        assert( internal.stats.isIntegerVals(Q,0,Inf) & length(Q) == R );

        % (3) If R = 0, then return an empty cell array.
        if R == 0
            if ~isempty(covariancepattern)
                warning(message('stats:LinearMixedModel:IgnoreCovariancePattern'));
            end
            covariancepattern = {};
            return;
        end           
        
        % (4) If covariancepattern is a single string or a single logical
        % matrix or a single double matrix, convert it into a cell array of
        % length 1.
        singlestring = internal.stats.isString(covariancepattern);
        singlelogical = islogical(covariancepattern) & ismatrix(covariancepattern);
        singledouble = isnumeric(covariancepattern) & ismatrix(covariancepattern);
        if singlestring || singlelogical || singledouble
            covariancepattern = {covariancepattern};           
        end
        
        % (5) Now covariancepattern must be a cell array of length R.
        if R > 1
            % <entry key="BadCovariancePattern_vector">''CovariancePattern'' must be a cell array of length {0}.</entry>            
            assertThat(     iscell(covariancepattern),'stats:LinearMixedModel:BadCovariancePattern_vector',num2str(R));
            assertThat(length(covariancepattern) == R,'stats:LinearMixedModel:BadCovariancePattern_vector',num2str(R));
        else
            % <entry key="BadCovariancePattern_scalar">''CovariancePattern'' must be a string, a logical matrix, or a cell array of length {0}.</entry>            
            assertThat(     iscell(covariancepattern),'stats:LinearMixedModel:BadCovariancePattern_scalar',num2str(R));
            assertThat(length(covariancepattern) == R,'stats:LinearMixedModel:BadCovariancePattern_scalar',num2str(R));
        end
       
        
        % (6) Loop over R and validate covariancepattern{i}.
        for i = 1:R            
            doubleinput = isnumeric(covariancepattern{i}) & ismatrix(covariancepattern{i});             
            % Convert double input to logical.
            if doubleinput
                covariancepattern{i} = logical(covariancepattern{i});
            end
             
             stringinput = internal.stats.isString(covariancepattern{i});
            logicalinput = islogical(covariancepattern{i}) & ismatrix(covariancepattern{i});             
             
            if stringinput
                if Q(i) == 0
                    covariancepattern{i} = internal.stats.getParamVal('FullCholesky',classreg.regr.LinearLikeMixedModel.AllowedCovariancePatterns,'CovariancePattern');
                else
                    covariancepattern{i} = internal.stats.getParamVal(covariancepattern{i},classreg.regr.LinearLikeMixedModel.AllowedCovariancePatterns,'CovariancePattern');
                end
            elseif logicalinput
                % <entry key="BadCovariancePatternElement_logical">Element {0} of ''CovariancePattern'' must be a square logical matrix of size {1}.</entry>               
                assertThat(all(size(covariancepattern{i}) == [Q(i),Q(i)]),'stats:LinearMixedModel:BadCovariancePatternElement_logical',num2str(i),num2str(Q(i)));                
            else
                % <entry key="BadCovariancePatternElement">Element {0} of ''CovariancePattern'' must be a string or a logical matrix.</entry>                
                assertThat(false,'stats:LinearMixedModel:BadCovariancePatternElement',num2str(i));
            end
        end
        
        % (7) Ensure size(covariancepattern,1) >= 1.
        if size(covariancepattern,1) == 1
            covariancepattern = covariancepattern';
        end
        
    end % end of validateCovariancePattern.
        
    function exclude = validateExclude(exclude,N)
%validateExclude - Validate the 'Exclude' parameter.   
%   exclude = validateExclude(exclude,N) takes exclude, a potential
%   'Exclude' parameter and an integer N and returns the validated exclude
%   parameter as a column vector or throws an error.
%
%   What is checked?
%
%   Either
%   (1) exclude is a vector containing integer values in the range 1 to N.
%
%       or
%
%   (2) exclude is a logical vector of length N.
%
%   If both (1) and (2) fail, then an error is thrown.

        % (1) N must be >= 0 and it must be an integer.
        assert(N >= 0 & internal.stats.isScalarInt(N));                 

        % (2) exclude as an integer vector with values from 1 to N.
        isintegervec = internal.stats.isIntegerVals(exclude,1,N); 
        
        % (3) exclude as a logical vector of length N.
        islogicalvec = isvector(exclude) & islogical(exclude) ...
            & length(exclude) == N;
        
        % (4) One of the following must be true. If not, throw error.
        % <entry key="BadExclude">The ''Exclude'' parameter must be a logical vector of length {0} or an integer vector with values in the range 1 to {1}.</entry>       
        assertThat(isintegervec | islogicalvec,'stats:LinearMixedModel:BadExclude',num2str(N),num2str(N));
        
        % (5) Ensure that exclude is a column vector.
        if size(exclude,1) == 1
            exclude = exclude';
        end
        
    end % end of validateExclude.
    
    function dummyvarcoding = validateDummyVarCoding(dummyvarcoding)
%validateDummyVarCoding - Validate the 'DummyVarCoding' parameter.    
%   dummyvarcoding = validateDummyVarCoding(dummyvarcoding) takes 
%   a tentative 'DummyVarCoding' parameter and returns the validated 
%   dummyvarcoding parameter or throws an error message.
%
%   What is checked?
%
%   (1) dummyvarcoding is a string that occurs in the cell array 
%   LinearMixedModel.AllowedDummyVarCodings.
%
%   If (1) is *not* satisfied, then an error message is thrown.

        dummyvarcoding = internal.stats.getParamVal(dummyvarcoding,...
            classreg.regr.LinearLikeMixedModel.AllowedDummyVarCodings,'DummyVarCoding');
    
    end % end of validateDummyVarCoding.

    function alpha = validateAlpha(alpha)
%validateAlpha - Validate the confidence level alpha.
%   alpha = validateAlpha(alpha) accepts a potential alpha value and
%   returns the validated alpha. If not valid, an error message is thrown.
%
%   What is checked?
%
%   (1) alpha must be a numeric, real, scalar.
%   (2) alpha must be >= 0 and <= 1.
% 
%   If (1) or (2) is not valid, an error message is thrown.
               
        % (1) Ensure that 'Alpha' is sensible.
        % <entry key="BadAlpha">The ''Alpha'' parameter must be a numeric, real, scalar between 0 and 1.</entry>
        assertThat(isnumeric(alpha) & isreal(alpha) & isscalar(alpha),'stats:LinearMixedModel:BadAlpha');               
        assertThat(        alpha >= 0 & alpha <= 1                   ,'stats:LinearMixedModel:BadAlpha');
       
    end % end of validateAlpha.

    function tf = validateLogicalScalar(tf,msgID)
%validateLogicalScalar - Validate a logical scalar input.
%   tf = validateLogicalScalar(tf,msgID) takes a potential logical scalar
%   value tf and validates it. If not valid, the error message with message
%   ID msgID is displayed. Output tf is always a logical scalar.
%
%   What is checked?
%
%   (1) tf must be a scalar.
%   (2) tf must be either true or 1.
%   (3) tf must be either false or 0.

        % (1) Ensure that msg is a string.
        assert( internal.stats.isString(msgID) );
        
        % (2) Now validate tf.
        if isscalar(tf)            
            if isnumeric(tf)
                if (tf == 1)
                    tf = true;
                elseif (tf == 0)
                    tf = false;
                end
            end                                                            
            assertThat(islogical(tf),msgID);
        else
           error(message(msgID)); 
        end 

    end % end of validateLogicalScalar.

    function wantconditional = validateConditional(wantconditional)
%validateConditional - Validate the parameter 'Conditional'.
%   wantconditional = validateConditional(wantconditional) takes a
%   potential value wantconditional for the 'Conditional' flag and
%   validates it. If not valid, an error message is thrown.
%
%   What is checked?
%
%   (1) wantconditional must be a scalar logical (true or false).
%
%   If (1) is not satisfied, then an error message is thrown.
        
        % <entry key="BadConditional">The ''Conditional'' parameter must be true (or 1) or false (or 0).</entry>        
        wantconditional = ...
            classreg.regr.LinearLikeMixedModel.validateLogicalScalar(wantconditional,'stats:LinearMixedModel:BadConditional');        
        
    end % end of validateConditional.
   
    function wantsimultaneous = validateSimultaneous(wantsimultaneous)
%validateSimultaneous - Validate the parameter 'Simultaneous'.
%   wantsimultaneous = validateSimultaneous(wantsimultaneous) takes a
%   potential value wantsimultaneous for the 'Simultaneous' flag and
%   validates it. If not valid, an error message is thrown.
%
%   What is checked?
%
%   (1) wantsimultaneous must be a scalar logical (true or false).
%
%   If (1) is not satisfied, then an error message is thrown.

        % <entry key="BadSimultaneous">The ''Simultaneous'' parameter must be true (or 1) or false (or 0).</entry>        
        wantsimultaneous = ...
            classreg.regr.LinearLikeMixedModel.validateLogicalScalar(wantsimultaneous,'stats:LinearMixedModel:BadSimultaneous');
        
    end % end of validateSimultaneous.
    
    function designtype = validateDesignType(designtype)
%validateDesignType - Validates the designtype input.
%   designtype = validateDesignType(designtype) accepts a potential value
%   of design type, validates it and returns it. If not valid, an error
%   message is thrown.
%
%   What is checked?
%
%   (1) designtype must be a string and must be either 'Fixed' or 'Random'.
%
%   If (1) is not satisfied, then an error is thrown.
     
        designtype = internal.stats.getParamVal(designtype,...
            {'Fixed','Random'},'DESIGNTYPE');
        
    end % end of validateDesignType.            

    function H = validateFEContrast(H,p)
%validateFEContrast - Validate the 'FEContrast' option in coefTest.        
%   H = validateFEContrast(H,p) takes a potential M-by-p contrast matrix H,
%   an integer p and returns the validated H. If not valid, an error
%   message is thrown.
%
%   What is checked?
%
%   (1) H is a numeric, real matrix with p columns.
%
%   If (1) does not hold, an error message is thrown.

        % (1) Ensure that p >= 0 and p is an integer.
        assert(p >= 0 & internal.stats.isScalarInt(p));
               
        % (2) Check H.
        % <entry key="BadFEContrast">H must be a numeric, real matrix with {0} columns.</entry>               
        assertThat(isnumeric(H) & isreal(H) & ismatrix(H) & size(H,2)==p,'stats:LinearMixedModel:BadFEContrast',num2str(p));

    end % end of validateFEContrast.
    
    function c = validateTestValue(c,M)
%validateTestValue - Validate the 'TestValue' option in coefTest.      
%   c = validateTestValue(c,M) takes a potential M-by-1 or 1-by-M vector c,
%   an integer M and returns the validated column vector c. If not valid, 
%   an error message is thrown.
%
%   What is checked?
%
%   (1) c is a numeric real vector with M elements.
%
%   If (1) does not hold, an error message is thrown.

        % (1) Ensure that M >= 0 and M is an integer.
        assert(M >= 0 & internal.stats.isScalarInt(M));

        % (2) Check c.
        % <entry key="BadTestValue">C must be a numeric, real vector of length {0}.</entry>       
        assertThat(isnumeric(c) & isreal(c) & isvector(c),'stats:LinearMixedModel:BadTestValue',num2str(M));
        if size(c,1) == 1
             c = c'; % Make column vector.
        end
        assertThat(       all( size(c)==[M,1] )          ,'stats:LinearMixedModel:BadTestValue',num2str(M));
        
    end % end of validateTestValue.
    
    function K = validateREContrast(K,M,q)
%validateREContrast - Validate the 'REContrast' option in coefTest.
%   K = validateREContrast(K,M,q) takes a potential M-by-q matrix K,
%   integers M and q and returns the validated K. If not valid, an error
%   message is thrown.
%
%   What is checked?
%
%   (1) K is a numeric, real matrix of size M-by-q.
%
%   If (1) does not hold, an error message is thrown.

        % (1) Ensure that M and q are >=0 and are integers.
        assert(M >= 0 & internal.stats.isScalarInt(M));
        assert(q >= 0 & internal.stats.isScalarInt(q));

        % (2) Check K.
        % <entry key="BadREContrast">''REContrast'' must be a numeric, real matrix of size {0} by {1}.</entry>              
        assertThat(isnumeric(K) & isreal(K) & ismatrix(K),'stats:LinearMixedModel:BadREContrast',num2str(M),num2str(q));   
        assertThat(          all(size(K)==[M,q])         ,'stats:LinearMixedModel:BadREContrast',num2str(M),num2str(q));   

    end % end of validateREContrast.
    
    function gnumbers = validateGNumbers(gnumbers,R)
%validateGNumbers - Validate the value of GNUMBERS positional input.      
%   gnumbers = validateGNumbers(gnumbers,R) accepts a potential integer
%   array gnumbers with elements between 1 and R and returns the validated
%   array gnumbers as a column vector.
%
%   What is checked?
%
%   (1) gnumbers is a numeric, real vector with elements between 1 and R.
%
%   If (1) is not satisfied, an error message is thrown.

        % (1) Ensure that R >= 0 and R is an integer.
        assert(R >= 0 & internal.stats.isScalarInt(R));

        % (2) Do the check.
        % <entry key="BadGNumbers">GNUMBERS must be an integer vector with elements between 1 and {0}.</entry>                
        assertThat(isnumeric(gnumbers) & isreal(gnumbers) & isvector(gnumbers),'stats:LinearMixedModel:BadGNumbers',num2str(R)); 
        assertThat(       internal.stats.isIntegerVals(gnumbers,1,R)          ,'stats:LinearMixedModel:BadGNumbers',num2str(R)); 

    end % end of validateGNumbers.

    function X = validateMatrix(X,XName,N,P)
%validateMatrix - Validates a numeric real matrix.       
%   X = validateMatrix(X,XName,N,P) accepts a potential matrix X of size 
%   N-by-p with name XName (a string) and validates it. If not valid, an 
%   error message is thrown.
%
%   What is checked?
%
%   (1) X is a numeric, real matrix.
%   (2) If N is given, we check that X has N rows.
%   (3) If P is given, we check that X has P columns.
%

        % (1) Set defaults for N and P.
        switch nargin        
            case 1
                % XName, N and P not given.
                XName = 'X';
                N = [];
                P = [];
            case 2
                % N, P not given.
                N = [];
                P = [];
            case 3
                % P not given.
                P = [];
        end
        
        % (2) Ensure that N, P >= 0 and N, P are integers.
        if ~isempty(N)
            assert(N >= 0 & internal.stats.isScalarInt(N));
        end
        if ~isempty(P)
            assert(P >= 0 & internal.stats.isScalarInt(P));
        end
        
        % (3) Ensure that XName is a string.
        assert( internal.stats.isString(XName) ); 
                
        % (4) Check if X is a numeric, real matrix.
        % <entry key="MustBeMatrix">{0} must be a numeric, real matrix.</entry>
        if islogical(X) && ismatrix(X)
            X = double(X);
        end
        assertThat(isnumeric(X) & isreal(X) & ismatrix(X),'stats:LinearMixedModel:MustBeMatrix',XName);        
        
        % (5) Check if X has N rows.
        if ~isempty(N)
            % <entry key="MustHaveRows">{0} must have {1} rows.</entry>            
            assertThat(size(X,1) == N,'stats:LinearMixedModel:MustHaveRows',XName,num2str(N));            
        end
        
        % (6) Check if X has P columns.
        if ~isempty(P)
            % <entry key="MustHaveCols">{0} must have {1} columns.</entry>           
            assertThat(size(X,2) == P, 'stats:LinearMixedModel:MustHaveCols',XName,num2str(P));            
        end    

    end % end of validateMatrix.
    
    function G = validateGroupingVar(G,GName,N)
%validateGroupingVar - Validate a grouping variable.      
%   G = validateGroupingVar(G,GName,N) accepts a potential grouping 
%   variable G with name GName expected to be of length N and returns the 
%   validated G as a column vector.       
%
%   What is checked?
%
%   (1) Either G is categorical, a cell array of strings, logical, or a 
%   character array. 
%       or
%   (2) G is a numeric, real vector.
%
%   (3) Is N is supplied, G must have length N.

        % (1) Set default values of GName and N.
        switch nargin
            case 1
                % GName and N not supplied.
                GName = 'G';
                N = [];
            case 2
                % N not supplied.
                N = [];
        end
        
        % (2) Ensure that N >= 0 and N is an integer.
        if ~isempty(N)
            assert(N >= 0 & internal.stats.isScalarInt(N));
        end
        
        % (3) Ensure that GName is a string.
        assert( internal.stats.isString(GName) );
        
        % (4) Check G.
        isdiscreteG = internal.stats.isDiscreteVar(G);
        isnumericG = isnumeric(G) & isreal(G) & isvector(G);
        
        % (5) Create error message.
        % <entry key="BadGroupingVar">{0} must be categorical, a cell array of strings, logical vector, a character array, or a numeric real vector.</entry>       
        assertThat(isdiscreteG | isnumericG,'stats:LinearMixedModel:BadGroupingVar',GName);
        
        % (6) Ensure G is a column vector.
        if size(G,1) == 1 && ~ischar(G)
            G = G';
        end
        
        % (7) Check size of G.
        if ~isempty(N)
            % <entry key="BadLength">{0} must have length {1}.</entry>                  
            assertThat(size(G,1) == N ,'stats:LinearMixedModel:BadLength',GName,num2str(N));
        end
        
    end % end of validateGroupingVar.
    
    function C = validateCellVectorOfStrings(C,CName,P,mustbeunique)
%validateCellVectorOfStrings - Validates a cell vector of strings.        
%   C = validateCellVectorOfStrings(C,CName,P,mustbeunique) takes a cell
%   vector of strings C of length P and a name CName for C and validate it.
%   C must have unique strings if mustbeunique is true, otherwise
%   repetitions are allowed. Default for mustbeunique is false. If not
%   valid, an error message is printed. The returned C is a column cell
%   array.
%
%   What is checked?
%
%   (1) C is a cell array of strings.
%   (2) C is a cell vector of length P.
%   (3) If mustbeunique is true, then C must have no repeated strings.

        % (1) Set default values for CName and P.
        switch nargin
            case 1
                % CName, P and mustbeunique not given.
                CName = 'C';
                P = [];
                mustbeunique = false;
            case 2
                % P and mustbeunique not given.
                P = [];
                mustbeunique = false;
            case 3
                % mustbeunique not given.
                mustbeunique = false;
        end        
        
        % (2) Ensure that C is a column cell vector of length P.
        C = classreg.regr.LinearLikeMixedModel.validateCellVector(C,CName,P);                
        
        % (3) Ensure that mustbeunique is true or false.
        assert( isscalar(mustbeunique) & islogical(mustbeunique) );
        
        % (4) Make sure C is a cell vector of strings.
        % <entry key="MustBeCellVectorOfStrings">{0} must be a cell vector of strings.</entry>        
        assertThat(iscellstr(C),'stats:LinearMixedModel:MustBeCellVectorOfStrings',CName);                                
        
        % (5) If mustbeunique is true, C must have unique strings.
        if mustbeunique == true
             n1 = length(C);
             n2 = length(unique(C));
            % <entry key="MustBeCellVectorOfStrings_unique">{0} must have unique strings.</entry>            
            assertThat(n1 == n2,'stats:LinearMixedModel:MustBeCellVectorOfStrings_unique',CName);            
        end
        
    end % end of validateCellVectorOfStrings.

    function S = validateString(S,SName)
%validateString - Validate a string.
%   S = validateString(S,SName) takes a potential string S with name SName
%   and returns the validated string. 
%
%   What is checked?
%
%   (1) S is a string.

        % (1) Ensure that SName is a string.
        assert( internal.stats.isString(SName) );
        
        % (2) Ensure that S is a string.
        % <entry key="MustBeString">{0} must be a string.</entry>        
        assertThat(internal.stats.isString(S),'stats:LinearMixedModel:MustBeString',SName);
        
    end % end of validateString.
    
    function C = validateCellVector(C,CName,R)
%validateCellVector - Validate that input is a cell vector.
%   C = validateCellVector(C,CName,R) takes a potential cell array C of
%   name CName and ensures that C is a cell vector of length R. The
%   returned C is a column cell vector.
%
%   What is checked?
%
%   (1) C is a cell vector.
%   (2) If R is not empty, C must have length R.

        % (1) Default values of CName and R.
        switch nargin
            case 1
                % CName and R not given.
                CName = 'C';
                R = [];
            case 2
                % R not given.
                R = [];
        end
            
        % (1) Ensure that CName is a string.
        assert( internal.stats.isString(CName) );
        
        % (2) Ensure that R >= 0 and R is an integer.
        if ~isempty(R)
            assert(R >= 0 & internal.stats.isScalarInt(R));
        end
        
        % (3) Ensure that C is a cell vector.
        % <entry key="MustBeCellVector">{0} must be a cell vector.</entry>        
        assertThat(iscell(C) & isvector(C),'stats:LinearMixedModel:MustBeCellVector',CName);
        
        % (4) Make C into a column cell vector.
        if size(C,1) == 1
            C = C';
        end
        
        % (5) Make sure C has length R.
        if ~isempty(R)
           % <entry key="MustBeCellVector_length">{0} must be a cell vector of length {1}.</entry>           
           assertThat(length(C) == R,'stats:LinearMixedModel:MustBeCellVector_length',CName,num2str(R));
        end        
        
    end % end of validateCellVector.  
    
    function ds = validateDataset(ds,dsName,dsref)
%validateDataset - Validate a dataset/table against a reference dataset/table.        
%   ds = validateDataset(ds,dsName,dsref) accepts a potential dataset/table ds, 
%   a string dsName and a reference dataset/table dsref and returns the validated 
%   ds.       
%
%   What is checked?
%
%   (1) ds must be a dataset/table.
%
%   (2) ds must have all variables that occur in dsref.
%
%   (3) The class of matching variables in ds and dsref must be the same.
%
%   If (1) or (2) or (3) is not satisfied, an error message is thrown.

        % (1) Ensure that dsName is a string.
        assert( internal.stats.isString(dsName) );
        
        if isa(dsref,'dataset') 
            dsref = dataset2table(dsref);
        end
        if isa(ds,'dataset')
            ds = dataset2table(ds);
        end
        % (2) Ensure that dsref is a dataset/table.
        assert( isa(dsref,'table') );

        % (3) Ensure that ds is a dataset/table.
        % <entry key="MustBeTable">{0} must be a dataset/table.</entry>        
        assertThat(isa(ds,'table'),'stats:LinearMixedModel:MustBeDataset',dsName);        
        
        % (4) Compare variable names in ds and dsref.
        varnames_dsref = dsref.Properties.VariableNames; 
        varnames_ds    =    ds.Properties.VariableNames;
        [tf,idx] = ismember( varnames_dsref, varnames_ds );                
        
        % (5) Make sure all variables in dsref are found in ds.
        % <entry key="Dataset_missingvarnames">Missing variable names in {0}.</entry>       
        assertThat(all(tf),'stats:LinearMixedModel:Dataset_missingvarnames',dsName);
        
        % (6) Compare classes of variables that matched.
        for j = 1:length(idx)    
            k = idx(j);
            if k ~= 0
                % varnames_dsref{j} is the same as varnames_ds{k}.
                   class_ds = class(ds.(varnames_ds{k}));
                class_dsref = class(dsref.(varnames_dsref{j}));
                % <entry key="Dataset_incorrectclass">Incorrect class for variable {0} in {1}.</entry>                
                assertThat(isequal(class_ds,class_dsref),'stats:LinearMixedModel:Dataset_incorrectclass',varnames_ds{k},dsName);
            end            
        end                

    end % end of validateDataset.
    
    function obj = validateObjectClass(obj,objName,className)
%validateObjectClass - Validates the class of an object.
%   obj = validateObjectClass(obj,objName,className) takes an object obj of
%   name given in string objName with expected class given in the string
%   className. The validated object obj is returned or an error is thrown
%   if obj does not have the class specified in className.
%
%   What is checked?
%
%   (1) isa(obj,className) should be true, otherwise we throw an error.

        % (1) Ensure that objName and className are strings.
        assert( internal.stats.isString(objName  ) );
        assert( internal.stats.isString(className) );
        
        % (2) Prepare the error message.
        % <entry key="BadClass">Object {0} should be of class {1}.</entry>       
        assertThat(isa(obj,className),'stats:LinearMixedModel:BadClass',objName,className);
        
    end % end of validateObjectClass.
    
    function checknesting = validateCheckNesting(checknesting)
%validateCheckNesting - Validate the parameter 'CheckNesting'.
%   checknesting = validateCheckNesting(checknesting) takes a potential 
%   value checknesting for the 'CheckNesting' flag and validates it. If not
%   valid, an error message is thrown.
%
%   What is checked?
%
%   (1) checknesting must be a scalar logical (true or false).
%
%   If (1) is not satisfied, then an error message is thrown.

        % <entry key="BadCheckNesting">The ''CheckNesting'' parameter must be true (or 1) or false (or 0).</entry>        
        checknesting = ...
            classreg.regr.LinearLikeMixedModel.validateLogicalScalar(checknesting,'stats:LinearMixedModel:BadCheckNesting');
        
    end % end of validateCheckNesting.
    
end

% Abstract static utility methods.
methods(Static,Abstract,Access='protected')
    
    fitmethod = validateFitMethod(fitmethod);
    w = validateWeights(w,N);
    [optimizer,optimizeroptions] = validateOptimizerAndOptions(optimizer,optimizeroptions);
    startmethod = validateStartMethod(startmethod);
    dfmethod = validateDFMethod(dfmethod);
    residualtype = validateResidualType(residualtype);
    verbose = validateVerbose(verbose);
    checkhessian = validateCheckHessian(checkhessian);
    
    checkNestingRequirement(smallModel,bigModel,smallModelName,bigModelName,isSimulatedTest);
    
end

% Redefined shared public methods from ParametricRegression.
methods (Access='public')
   
    function [feci,reci] = coefCI(model,varargin)
%coefCI Confidence intervals for coefficients.
%   FECI = coefCI(LME) computes 95% confidence intervals for the fixed
%   effects parameters in the linear mixed effects model LME. The output
%   FECI is a P-by-2 matrix where P is the number of fixed effects
%   parameters in LME. Rows of FECI from top to bottom correspond
%   respectively to the P-by-1 fixed effects vector BETA displayed from top
%   to bottom in the tabular display from the fixedEffects method. Column 1
%   of FECI displays lower confidence limits and column 2 of FECI displays
%   upper confidence limits.
%
%   [FECI,RECI] = coefCI(LME) also returns 95% confidence intervals for
%   random effects parameters in LME. The output RECI is a Q-by-2 matrix
%   where Q is the total number of random effects parameters in LME. Rows
%   of RECI from top to bottom correspond respectively to the Q-by-1 random
%   effects vector B displayed from top to bottom in the tabular display
%   from the randomEffects method. Column 1 of RECI displays lower
%   confidence limits and column 2 of RECI displays upper confidence
%   limits.
%
%   [FECI,RECI] = coefCI(LME,'PARAM','VALUE',...) accepts optional
%   name/value pairs:
%      
%           'Name'     'Value'
%          'Alpha'     ALPHA, a number between 0 and 1. Computes 
%                      100*(1-ALPHA)% confidence intervals. Default is
%                      ALPHA=0.05 for 95% confidence intervals.
%  
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate degrees of freedom DF for confidence
%                      interval calculation. Options are 'Satterthwaite',
%                      'Residual' and 'None'. If 'DFMethod' is
%                      'Satterthwaite', a Satterthwaite approximation is
%                      used to compute DF. If 'DFMethod' is 'Residual', the
%                      DF values are assumed to be constant and equal to
%                      (N-P) where N is the number of observations and P is
%                      the number of fixed effects. If 'DFMethod' is
%                      'None', then all DF values are set to infinity.
%                      Default is 'Residual'.
%
%   Example: Fit a model with Weight as a fixed effect, and a random effect
%            due to Model_Year. Compute confidence intervals for the
%            intercept and the coefficient of Weight.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
%      coefCI(lme)
%
%   See also coefTest, fixedEffects, randomEffects, covarianceParameters.

        % (1) Default parameter values.
        dfltDFMethod = 'Residual';
           dfltAlpha = 0.05;

        % (2) Optional parameter names and their default values.
        paramNames =   {'DFMethod',   'Alpha'};
        paramDflts = {dfltDFMethod, dfltAlpha};
           
        % (3) Parse optional parameter name/value pairs.
        [dfmethod,alpha] ...
            = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % (4) Validate optional parameter values.
        dfmethod = model.validateDFMethod(dfmethod);
           alpha = model.validateAlpha(alpha);

        % (5) First compute feci.
        fetable = fixedEffects(model.slme,alpha,dfmethod);
        feci    = [fetable.Lower, fetable.Upper];
                    
        % (6) Compute reci if asked.
        if nargout > 1
            
            retable = randomEffects(model.slme,alpha,dfmethod);
            reci    = [retable.Lower, retable.Upper];
            
        end

    end % end of coefCI.
    
%     function [feci,reci] = coefCIOLD(model,varargin)
% %coefCI Confidence intervals for coefficients.
% %   FECI = coefCI(LME) computes 95% confidence intervals for the fixed
% %   effects parameters in the linear mixed effects model LME. The output
% %   FECI is a P-by-2 matrix where P is the number of fixed effects
% %   parameters in LME. Rows of FECI from top to bottom correspond
% %   respectively to the P-by-1 fixed effects vector BETA displayed from top
% %   to bottom in the tabular display from the fixedEffects method. Column 1
% %   of FECI displays lower confidence limits and column 2 of FECI displays
% %   upper confidence limits.
% %
% %   [FECI,RECI] = coefCI(LME) also returns 95% confidence intervals for
% %   random effects parameters in LME. The output RECI is a Q-by-2 matrix
% %   where Q is the total number of random effects parameters in LME. Rows
% %   of RECI from top to bottom correspond respectively to the Q-by-1 random
% %   effects vector B displayed from top to bottom in the tabular display
% %   from the randomEffects method. Column 1 of RECI displays lower
% %   confidence limits and column 2 of RECI displays upper confidence
% %   limits.
% %
% %   [FECI,RECI] = coefCI(LME,'PARAM','VALUE',...) accepts optional
% %   name/value pairs:
% %      
% %           'Name'     'Value'
% %          'Alpha'     ALPHA, a number between 0 and 1. Computes 
% %                      100*(1-ALPHA)% confidence intervals. Default is
% %                      ALPHA=0.05 for 95% confidence intervals.
% %  
% %       'DFMethod'     Specifies the method to use for computing the
% %                      approximate degrees of freedom DF for confidence
% %                      interval calculation. Options are 'Satterthwaite',
% %                      'Residual' and 'None'. If 'DFMethod' is
% %                      'Satterthwaite', a Satterthwaite approximation is
% %                      used to compute DF. If 'DFMethod' is 'Residual', the
% %                      DF values are assumed to be constant and equal to
% %                      (N-P) where N is the number of observations and P is
% %                      the number of fixed effects. If 'DFMethod' is
% %                      'None', then all DF values are set to infinity.
% %                      Default is 'Residual'.
% %
% %   Example: Fit a model with Weight as a fixed effect, and a random effect
% %            due to Model_Year. Compute confidence intervals for the
% %            intercept and the coefficient of Weight.
% %      load carsmall
% %      ds = dataset(MPG,Weight,Model_Year);
% %      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
% %      coefCI(lme)
% %
% %   See also coefTest, fixedEffects, randomEffects, covarianceParameters.
% 
%         % (1) Default parameter values.
%         dfltDFMethod = 'Residual';
%            dfltAlpha = 0.05;
% 
%         % (2) Optional parameter names and their default values.
%         paramNames =   {'DFMethod',   'Alpha'};
%         paramDflts = {dfltDFMethod, dfltAlpha};
%            
%         % (3) Parse optional parameter name/value pairs.
%         [dfmethod,alpha] ...
%             = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   
% 
%         % (4) Validate optional parameter values.
%         dfmethod = model.validateDFMethod(dfmethod);
%            alpha = model.validateAlpha(alpha);
% 
%         % (5) First compute feci.           
%             % (1) How many fixed effects in the model?
%             p = model.slme.p;
%             
%             % (2) Loop 1:p and call betaCI in StandardLinearMixedModel. If
%             % beta is the fixed effects vector, betaCI gets CIs on c'*beta.
%             % CIs on beta(i) can be obtained by using the contrast vector c
%             % such that:
%             %           c = zeros(p,1);
%             %           c(i) = 1;
%             % i.e., there is a 1 in the i th position and 0's elsewhere.
%             feci = zeros(p,2);
%             for i = 1:p
%                 c = zeros(p,1);
%                 c(i) = 1;
%                 feci(i,:) = betaCI(model.slme,c,alpha,dfmethod); 
%             end
%             
%         % (6) Compute reci if asked.
%         if nargout > 1
%            
%             % (1) Store the joint covariance of [betaHat - beta, bHat - b]
%             % in the slme object.
%             model.slme = storeCovBetaHatBHat(model.slme);            
%             
%             % (2) How many random effects in the model?
%             q = model.slme.q;
%             
%             % (3) Loop 1:q and call betaBCI in StandardLinearMixedModel. If
%             % beta is the fixed effects vector and b is the random effects 
%             % vector, betaBCI gets CIs on t'*beta + s'*b. For CIs on b(i),
%             % set t = zeros(p,1) and create a vector s such that:
%             %           s = zeros(q,1);
%             %           s(i) = 1;
%             % i.e., there is a 1 in the i th position and 0's elsewhere.            
%             reci = zeros(q,2);
%             t = zeros(p,1);
%             for i = 1:q
%                 s = zeros(q,1);
%                 s(i) = 1;
%                 reci(i,:) = betaBCI(model.slme,t,s,alpha,dfmethod);
%             end
%                         
%             % (4) Unstore the covariance of [betaHat - beta, bHat - b] from
%             % the object.
%             model.slme = unstoreCovBetaHatBHat(model.slme); 
%             
%         end
% 
%     end % end of coefCIOLD.

    function [P,F,DF1,DF2] = coefTest(model,H,c,varargin)
%coefTest Linear hypothesis test on coefficients.
%   PVAL = coefTest(LME) computes the p-value for an F test that all fixed
%   effects parameters in the linear mixed effects model LME except the
%   intercept are zero.
%
%   PVAL = coefTest(LME,H) computes the p-value for an F test on the fixed
%   effects part of LME using the M-by-P matrix H where P is the number of
%   fixed effects parameters in LME. Each row of H represents 1 contrast
%   and the columns of H from left to right correspond respectively to the
%   P-by-1 fixed effects vector BETA displayed from top to bottom in the
%   tabular display from the fixedEffects method. The output PVAL is the
%   p-value for an F test that H*BETA = 0. To include contrasts that
%   involve the random effects, use the 'REContrast' parameter.
%
%   PVAL = coefTest(LME,H,C) also specifies a M-by-1 vector C for testing 
%   the hypothesis H*BETA = C.
%
%   PVAL = coefTest(LME,H,C,PARAM1,VALUE1,...) accepts optional name/value
%   pairs to control the calculation of PVAL.
%
%             Name     Value
%       'DFMethod'     A string that specifies the method to use for
%                      computing the approximate denominator degrees of
%                      freedom DF for the F test. If 'DFMethod' is
%                      'Satterthwaite', a Satterthwaite approximation is
%                      used to compute DF. If 'DFMethod' is 'Residual', the
%                      DF value is assumed to be equal to (N-P) where N is
%                      the number of observations and P is the number of
%                      fixed effects. If 'DFMethod' is 'None', then the DF
%                      value is taken to be infinity. Default is
%                      'Residual'.
%
%       'REContrast'   A M-by-Q matrix K where Q is the number of random 
%                      effects parameters in LME. Each row of K represents
%                      1 contrast and the columns of K from left to right
%                      correspond respectively to the Q-by-1 random effects
%                      vector B displayed from top to bottom in the tabular
%                      display from the randomEffects method. The output
%                      PVAL is the p-value for an F test H*BETA + K*B = C.
%
%   [PVAL,F,DF1,DF2] = coefTest(...) also returns the F-statistic F, the
%   numerator degrees of freedom DF1 for F, and the denominator degrees of
%   freedom DF2 for F. DF1 is equal to the number of linearly independent
%   rows in H, or [H,K] depending on the call to coefTest. The value of DF2
%   depends on the option selected for 'DFMethod'.
%
%   Example: Fit a model with two fixed-effect predictors and a random
%            effect. Test for the significance of the Cylinders term. The
%            p-value is the same as shown in the anova table.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year,Cylinders);
%      ds.Cylinders = nominal(ds.Cylinders);
%      lme = fitlme(ds,'MPG ~ Weight + Cylinders + (1|Model_Year)')
%      p = coefTest(lme,[0 0 1 0;0 0 0 1])
%      anova(lme)
%
%      % Test for no difference between 6-cylinder and 8-cylinder cars
%      coefTest(lme,[0 0 1 -1])
%
%   See also ANOVA, coefCI, fixedEffects, randomEffects, covarianceParameters. 

        % (1) Default parameter values.        
        dfltREContrast = [];
          dfltDFMethod = 'Residual';

        % (2) Optional parameter names and their default values.
        paramNames =   {'REContrast',   'DFMethod'};
        paramDflts = {dfltREContrast, dfltDFMethod};
           
        % (3) Parse optional parameter name/value pairs.
        [K,dfmethod] = ...
            internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % (4) Get the number of fixed and random coefficients.
        p = model.slme.p;
        q = model.slme.q;                
        
        if nargin < 2
            % (5) If H is not supplied, fill in default value. Default is
            % to test all fixed effects coefficients except the intercept
            % against zero. Initialize H to be the identity matrix to test
            % all fixed effects coefficients and then delete the row of H
            % corresponding to the intercept coefficient.
            H = eye(p);
            % Find the index of the intercept "term".
            interceptTermIdx = ...
                find(all(model.Formula.FELinearFormula.Terms == 0,2));
            % Which coefficient index is the intercept "term"?
            if isscalar(interceptTermIdx)
                interceptCoeffIdx = ...
                    model.FixedInfo.XCols2Terms == interceptTermIdx;
                % Delete rows in H corresponding to the intercept.
                H(interceptCoeffIdx,:) = [];
            end
        end
        
        if nargin < 3                               
            % (6) If c is not supplied, fill in default value. Default is
            % to test against zero.
            c = zeros(size(H,1),1);
        end
                
        % (7) Validate H, c and dfmethod.
        H = model.validateFEContrast(H,p);
        M = size(H,1);        
        c = model.validateTestValue(c,M);        
        dfmethod = model.validateDFMethod(dfmethod);

        % (8) If K is empty, we just want a test on fixed effects,
        % otherwise we need to validate K.
        if isempty(K)
            % Test H*beta = c via model.slme.
            [P,F,DF1,DF2] = betaFTest(model.slme,H,c,dfmethod);
        else
            % Validate K.
            K = model.validateREContrast(K,M,q);             
            % Test H*beta + K*b = c via model.slme.
            [P,F,DF1,DF2] = betaBFTest(model.slme,H,K,c,dfmethod);
        end

    end % end of coefTest.
    
end

% Shared public methods of LinearLikeMixedModel.
methods (Access='public')
    
    function [D,gnames] = designMatrix(model,designtype,gnumbers)
%designMatrix Extracts the fixed or random effects design matrices.
%   X = designMatrix(LME) or designMatrix(LME,'Fixed') returns the N-by-P
%   fixed effects design matrix X for the linear mixed effects model LME
%   where N is the number of observations and P is the number of fixed
%   effects.
%
%   D = designMatrix(LME,'Random') returns the overall random effects
%   design matrix corresponding to a vector B of all random effects in the
%   linear mixed effects model LME. Suppose LME has R grouping variables
%   named g_1,...,g_R. Let Q_1,...,Q_R be the length of random effects
%   vectors associated with g_1,...,g_R respectively. Also, suppose
%   g_1,...,g_R have levels M_1,...,M_R respectively. Then B will be a
%   column vector of length Q_1*M_1 + Q_2*M_2 + ... + Q_R*M_R. B is made by
%   concatenating the BLUPs of random effects vectors corresponding to each
%   level of each grouping variable in the following order:
%
%   g_1 level 1, g_1 level 2,..., g_1 level M_1 followed by
%   g_2 level 1, g_2 level 2,..., g_2 level M_2 followed by
%      .            .                .           
%   g_R level 1, g_R level 2,..., g_R level M_R
%
%   If LME has N observations then D will be of size N-by-length(B) such
%   that D*B is a N-by-1 vector that represents the contribution of all
%   random effects to the response of the LME.
%
%   DSUB = designMatrix(LME,'Random',GNUMBERS) returns a submatrix of the
%   full random effects design matrix. GNUMBERS is a length K integer array
%   with elements in the range 1 to R. DSUB is a subset of the full random
%   effects design matrix corresponding to the grouping variable names
%   indicated by integers in GNUMBERS. For example, suppose GNUMBERS is
%   [1,R] then this specifies only grouping variables g_1 and g_R. Let BSUB
%   be a vector made by concatenating BLUPs of random effects vectors
%   corresponding to each level of g_1 and g_R in the following order:
%
%   g_1 level 1, g_1 level 2,..., g_1 level M_1 followed by         
%   g_R level 1, g_R level 2,..., g_R level M_R
%
%   Then DSUB will be a N-by-length(BSUB) matrix such that DSUB*BSUB
%   represents the contribution of all random effects corresponding to
%   grouping variables g_1 and g_R to the response of the LME. If GNUMBERS
%   is empty, the full random effects design matrix is returned.
%
%   [DSUB,GNAMES] = designMatrix(LME,DESIGNTYPE,GNUMBERS) also returns a 
%   K-by-1 cell array containing the names of grouping variables
%   corresponding to integers in GNUMBERS if DESIGNTYPE is 'Random'. If
%   DESIGNTYPE is 'Fixed' then GNAMES is [] and GNUMBERS is ignored.
%        
%   Example: Fit a quadratic model with continuous and categorical fixed
%            effects. Look at a few rows of the design matrix showing the
%            constant term, the Weight term, and dummy variables for the
%            Origin term.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year,Origin);
%      lme = fitlme(ds,'MPG ~ Weight + Origin + (1|Model_Year)');
%      dm = designMatrix(lme);
%      disp(dm(25:30,:))
%
%   See also RESPONSE, FITTED, RESIDUALS.

        % (1) Fill in defaults for positional parameters if needed.
        narginchk(1,3);
        switch nargin
            case 1
                designtype = 'Fixed';
                gnumbers = [];
            case 2
                gnumbers = [];
        end
                
        % (2) Validate designtype.
        designtype = model.validateDesignType(designtype);                        
        
        % (3) Handle the 2 cases.
        switch lower(designtype)
            case 'fixed'
                D = fixedEffectsDesign(model);
                gnames = [];
            case 'random'
                [D,gnames] = randomEffectsDesign(model,gnumbers);                
        end
        
    end % end of designMatrix.
            
    function [beta,betanames,fetable] = fixedEffects(model,varargin)
%fixedEffects Returns estimates of fixed effects and related statistics.
%   BETA = fixedEffects(LME) returns a vector of estimated fixed effects
%   from a fitted linear mixed effects model LME.
%
%   [BETA,BETANAMES] = fixedEffects(LME) also returns a dataset array
%   BETANAMES containing the name of each fixed effects coefficient in
%   BETA.
%
%   [BETA,BETANAMES,STATS] = fixedEffects(LME) also returns a dataset array
%   STATS containing the estimates of fixed effects and related statistics.
%   Table STATS has one row for each fixed effect and the following
%   columns:
%
%       Name        Name of the fixed effect coefficient
%       Estimate    Estimated coefficient value
%       SE          Standard error of the estimate
%       tStat       t statistic for a test that the coefficient is zero
%       DF          Estimated degrees of freedom for the t statistic
%       pValue      p-value for the t statistic
%       Lower       Lower limit of a 95% confidence interval
%       Upper       Upper limit of a 95% confidence interval
%
%   [BETA,BETANAMES,STATS] = fixedEffects(LME,'PARAM','VALUE',...) also
%   specifies optional parameter name/value pairs to control the fixed
%   effects statistics:
%  
%           'Name'     'Value'
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate degrees of freedom (DF) for the t
%                      statistics that test the fixed effects coefficients
%                      against 0. Options are 'Satterthwaite','Residual'
%                      and 'None'. If 'DFMethod' is 'Satterthwaite', a
%                      Satterthwaite approximation is used to compute DF.
%                      If 'DFMethod' is 'Residual', the DF values are
%                      assumed to be constant and equal to (N-P) where N is
%                      the number of observations and P is the number of
%                      fixed effects. If 'DFMethod' is 'None', then all DF
%                      values are set to infinity. Default is 'Residual'.
%
%          'Alpha'     ALPHA, a number between 0 and 1. The columns 'Lower'
%                      and 'Upper' in table STATS will correspond to the
%                      lower and upper limits respectively of 100*(1-ALPHA)
%                      confidence intervals for fixed effects. Default is
%                      ALPHA=0.05 for 95% confidence intervals.
%
%   Example: Fit a model with a fixed effect and a random effect. Get the
%            estimated fixed effects (intercept and slope).
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
%      fixedEffects(lme)
%
%   See also coefCI, coefTest, randomEffects.

        % (1) Default parameter values.
        dfltDFMethod = 'Residual';
           dfltAlpha = 0.05;

        % (2) Optional parameter names and their default values.
        paramNames =   {'DFMethod',   'Alpha'};
        paramDflts = {dfltDFMethod, dfltAlpha};
           
        % (3) Parse optional parameter name/value pairs.
        [dfmethod,alpha] ...
            = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % (4) Validate optional parameter values.
        dfmethod = model.validateDFMethod(dfmethod);
           alpha = model.validateAlpha(alpha);
        
        % (5) Get beta.
        beta = model.slme.betaHat;
        
        % (6) Get betanames if required.
        if nargout > 1
            betanames = table(model.FixedInfo.XColNames','VariableNames',{'Name'});
        end
        
        % (7) Call fixedEffects methods on model.slme to get the stats
        % table if required.
        if nargout > 2            
            fetable = fixedEffects(model.slme,alpha,dfmethod);
            
            % fetable doesn't have fixed effects variable names. Add this
            % info as a new column with name: 'Name'.
            fetable = [betanames,fetable];
            
            % Put an informative title on fetable.
            % <entry key="Title_fetable">Fixed effect coefficients: DFMethod = ''{0}'', Alpha = {1}</entry>
            ttl = getString(message('stats:LinearMixedModel:Title_fetable',dfmethod,num2str(alpha)));
            fetable = classreg.regr.lmeutils.titleddataset(fetable,ttl);
        end

    end % end of fixedEffects.

    function [b,bnames,retable] = randomEffects(model,varargin)
%randomEffects Extract estimates of random effects and related statistics.
%   B = randomEffects(LME) returns estimates of the best linear unbiased
%   predictors (BLUPs) of all random effects in the linear mixed effects
%   model LME. Suppose LME has R grouping variables named g_1,...,g_R. Let
%   Q_1,...,Q_R be the length of random effects vectors associated with
%   g_1,...,g_R respectively. Also, suppose g_1,...,g_R have levels
%   M_1,...,M_R respectively. Then B will be a column vector of length
%   Q_1*M_1 + Q_2*M_2 + ... + Q_R*M_R. B is made by concatenating the BLUPs
%   of random effects vectors corresponding to each level of each grouping
%   variable in the following order:
%
%   g_1 level 1, g_1 level 2,..., g_1 level M_1 followed by
%   g_2 level 1, g_2 level 2,..., g_2 level M_2 followed by
%      .            .                .           
%   g_R level 1, g_R level 2,..., g_R level M_R
%
%   [B,BNAMES] = randomEffects(LME) also returns a table array BNAMES
%   containing the names of the coefficients in B.
%
%   [B,BNAMES,STATS] = randomEffects(LME) also returns a table array
%   STATS containing the estimates of random effects and related
%   statistics. Table STATS has one row for each random effect associated
%   with a particular grouping variable, level and predictor name with the
%   following columns:
%
%       Group       Grouping variable associated with this random effect
%       Level       Level within the grouping variable
%       Name        Name of the random effect coefficient
%       Estimate    BLUP of random effect
%       SEPred      Standard error of (BLUP minus the random effect)
%       tStat       t statistic for a test that the random effect is zero
%       DF          Estimated degrees of freedom for the t statistic
%       pValue      p-value for the t statistic
%       Lower       Lower limit of a 95% confidence interval
%       Upper       Upper limit of a 95% confidence interval
%
%   [B,BNAMES,STATS] = randomEffects(LME,'PARAM','VALUE',...) specifies
%   optional parameter name/value pairs to control the calculation of
%   random effects statistics:
%  
%           'Name'     'Value'
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate degrees of freedom (DF) for the t
%                      statistics that test the random effects coefficients
%                      against 0. Options are 'Satterthwaite','Residual'
%                      and 'None'. If 'DFMethod' is 'Satterthwaite', a
%                      Satterthwaite approximation is used to compute DF.
%                      If 'DFMethod' is 'Residual', the DF values are
%                      assumed to be constant and equal to (N-P) where N is
%                      the number of observations and P is the number of
%                      fixed effects. If 'DFMethod' is 'None', then all DF
%                      values are set to infinity. Default is 'Residual'.
%
%          'Alpha'     ALPHA, a number between 0 and 1. The columns 'Lower'
%                      and 'Upper' in table STATS will correspond to the
%                      lower and upper limits respectively of 100*(1-ALPHA)
%                      confidence intervals for random effects. Default is
%                      ALPHA=0.05 for 95% confidence intervals.
%
%   Example: Fit a model with a fixed effect and random effects. Get the
%            estimated random effects. They are highly negatively
%            correlated. Verify this by plotting them.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (Weight|Model_Year)')
%      re = randomEffects(lme)
%      plot(re(1:2:end),re(2:2:end),'rs')
%
%   See also coefCI, coefTest, fixedEffects.

        % (1) Default parameter values.
        dfltDFMethod = 'Residual';
           dfltAlpha = 0.05;

        % (2) Optional parameter names and their default values.
        paramNames =   {'DFMethod',   'Alpha'};
        paramDflts = {dfltDFMethod, dfltAlpha};
           
        % (3) Parse optional parameter name/value pairs.
        [dfmethod,alpha] ...
            = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % (4) Validate optional parameter values.
        dfmethod = model.validateDFMethod(dfmethod);
           alpha = model.validateAlpha(alpha);
        
        % (5) Get estimated BLUP of b.
        b = model.slme.bHat;
           
        % (6) Get names of elements of b in a table bnames.
        if nargout > 1
            bnames = model.RandomInfo.ZsColNames;
        end
        
        % (7) Call randomEffects methods on model.slme to get the stats
        % table if required.
        if nargout > 2            
            retable = randomEffects(model.slme,alpha,dfmethod);
            
            % retable doesn't have random effects name info. Add this
            % info by appending model.RandomInfo.ZsColNames in front of
            % retable.
            retable = [bnames,retable];
            
            % Put an informative title on retable.
            % <entry key="Title_retable">Random effect coefficients: DFMethod = ''{0}'', Alpha = {1}</entry>
            ttl = getString(message('stats:LinearMixedModel:Title_retable',dfmethod,num2str(alpha)));
            retable = classreg.regr.lmeutils.titleddataset(retable,ttl);            
        end

    end % end of randomEffects.

    function [PSI,mse,covtable] = covarianceParameters(model,varargin)
%covarianceParameters Extract estimated covariance parameters of the model.
%   PSI = covarianceParameters(LME) extracts the estimated covariance
%   parameters that parameterize the prior covariance of random effects. If
%   the linear mixed effects model LME has R grouping variables named
%   g_1,...,g_R then the output PSI will be a R-by-1 cell array such that
%   PSI{i} will contain the covariance matrix of random effects associated
%   with grouping variable g_i. The order in which grouping variables are
%   assigned numbers 1 to R is the same order in which grouping variables
%   are entered into the FITLME or FITLMEMATRIX functions.
%
%   [PSI,MSE] = covarianceParameters(LME) also extracts an estimate of the
%   residual variance.
%
%   [PSI,MSE,STATS] = covarianceParameters(LME) also returns a cell array
%   STATS of length (R+1) containing covariance parameters and related
%   statistics. STATS{i} is a table array containing statistics on
%   covariance parameters for the i-th grouping variable. STATS{R+1}
%   contains statistics on the residual standard deviation. STATS{i} 
%   contains columns that name each covariance parameter as well as the 
%   following columns:
%
%       Group       Grouping variable name
%       Estimate    Estimate of the covariance parameter
%       Lower       Lower limit of a 95% confidence interval
%       Upper       Upper limit of a 95% confidence interval
%
%   NOTE: It is recommended that the presence or absence of covariance
%   parameters in LME be tested using the COMPARE method, which uses a
%   likelihood ratio test.
%
%   [PSI,MSE,STATS] = covarianceParameters(LME,PARAM1,VALUE1,...)
%   specifies optional name/value pairs as follows:
%
%           'Name'     'Value'
%          'Alpha'     ALPHA, a number between 0 and 1. The columns 'Lower'
%                      and 'Upper' in table STATS{i} will correspond to
%                      the lower and upper limits respectively of
%                      100*(1-ALPHA) confidence intervals for covariance
%                      parameters. Default is ALPHA = 0.05 for 95% 
%                      confidence intervals. 
%
%   Example: Fit a model with two correlated random effects, and get their
%            estimated covariance matrix.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (Weight|Model_Year)')
%      V = covarianceParameters(lme);
%      V{1}
%
%   See also COMPARE, fixedEffects, randomEffects.

        % (1) Default parameter values.
        dfltAlpha = 0.05;

        % (2) Optional parameter names and their default values.
        paramNames = {  'Alpha', 'WantCIs'};
        paramDflts = {dfltAlpha, true};
           
        % (3) Parse optional parameter name/value pairs.
        [alpha,wantCIs] = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % (4) Validate optional parameter values.
        alpha = model.validateAlpha(alpha);
        assert( isscalar(wantCIs) & islogical(wantCIs) );

        % (5) How many grouping variables do we have?
        R = model.GroupingInfo.R;
        
        % (6) Get mse.
        mse = (model.slme.sigmaHat)^2;

        % (7) First compute PSI.
        PSI = cell(R,1);
        for k = 1:R           
            % (1) model.slme.Psi is always a BlockedCovariance object. Get
            % the lower triangular Cholesky factor from the covariance
            % matrix for the k th grouping variable.
            matk = model.slme.Psi.Matrices{k};
            Lk = getLowerTriangularCholeskyFactor(matk);  
            
            % (2) Now form the full covariance matrix for k th grouping
            % variable including mse.
            PSI{k} = mse*(Lk*Lk');            
        end
        
        % (8) Get covtable if required.
        if nargout > 2            
            % (1) Forward the call to model.slme. tbl is a table that 
            % always has a column Estimate but it may not have columns 
            % Lower and Upper if CI calculation cannot be performed. The
            % number of rows in tbl is equal to 1 + the total number of
            % canonical parameters (i.e., parameters of interest) for all
            % grouping variables. The last row of tbl corresponds to the 
            % residual error term. 
            tbl = covarianceParameters(model.slme,alpha,wantCIs);
            
            % (2) Get names of canonical parameters in all but the last row
            % of tbl.
            %
            % names should be a R-by-1 cell array of datasets if R >= 1. If
            %
            %   matk = model.slme.Psi.Matrices{k} then 
            %
            %   size(names{k},1) = matk.NumParametersExcludingSigma
            %
            % Also \sum_k{ size(names{k},1) } + 1 = size(tbl,1)
            if R >= 1
                names = getCanonicalParameterNames(model.slme.Psi);
            else
                names = cell(0,1);
            end
            
            % (3) Suppose the LME model has 2 grouping variables. The first
            % grouping variable has 2 canonical parameters and the second
            % grouping variable has 1 canonical parameter. In this case the
            % table tbl will look something like this:
            %
            %   Estimate    Lower   Upper
            %      x          x       x   <- names{1}(1,:)
            %      x          x       x   <- names{1}(2,:)
            %      x          x       x   <- names{2}(1,:)
            %      x          x       x   <- Res Std
            %
            % The markup on the right of Upper shows that names{1} contains
            % the name info for the first 2 rows of tbl and names{2}
            % contains info for the 3 rd row of tbl. The last row of tbl is
            % named 'Res Std' since it is an estimate of the residual
            % standard deviation.
            %
            % Build the output array covtable. covtable{k} should contain
            % the name info from names{k} and the relevant rows from tbl
            % for the k th grouping variable.        
            covtable = cell(R+1,1);
            offset = 0;
            for k = 1:R            
                % (1) Build an index idxk to extract relevant rows from tbl
                % for grouping variable k.
                  matk = model.slme.Psi.Matrices{k};
                startk = offset + 1;
                  endk = offset + matk.NumParametersExcludingSigma;
                  idxk = startk : endk;
                  
                % (2) k the element of covtable.  
                covtable{k} = [names{k},tbl(idxk,:)];
                
                % (3) Add a title for covtable{k}.
                %ttl = ['Group: ',model.GroupingInfo.GNames{k},...
                %    ', Covariance Type: ',model.slme.Psi.Matrices{k}.Type];
                % <entry key="Title_covtable">Covariance Type: {0}</entry>
                ttl = getString(message('stats:LinearMixedModel:Title_covtable',model.slme.Psi.Matrices{k}.Type));
                covtable{k} = classreg.regr.lmeutils.titleddataset(covtable{k},ttl);
                
                % (4) Update offset to go to the next grouping variable.
                offset = endk;
            end                         
            
            % (4) The (R+1) st row of covtable should contain info on the
            % residual standard deviation. There is no grouping variable
            % for the residual error term but let's indicate that by the
            % word 'Error'. 
            if isequal(class(model),'GeneralizedLinearMixedModel')
                res_std_name = 'sqrt(Dispersion)';
            else
                res_std_name = 'Res Std';
            end
            covtable{R+1} = [classreg.regr.lmeutils.titleddataset({'Error','Group'}, {{res_std_name},'Name'}),...
                classreg.regr.lmeutils.titleddataset(tbl(end,:))];  
            covtable{R+1}.Group = char(covtable{R+1}.Group);
        end

    end % end of covarianceParameters.    
    
    function hout = plotResiduals(model,plottype,varargin)
% plotResiduals Plot residuals of fitted model
%   plotResiduals(LME,PLOTTYPE) plots the raw conditional residuals from
%   model LME in a plot of type PLOTTYPE. Valid values for PLOTTYPE are:
%  
%      'caseorder'     residuals vs. case (row) order
%      'fitted'        residuals vs. fitted values
%      'histogram'     histogram (default)
%      'lagged'        residual vs. lagged residual (r(t) vs. r(t-1))
%      'probability'   normal probability plot
%      'symmetry'      symmetry plot
%  
%   plotResiduals(MODEL,PLOTTYPE,'PARAM','VALUE',...) accepts optional
%   parameter name/value pairs:
%
%           'Name'     'Value'
%    'ResidualType'     Valid values are 'Raw' (default), 'Standardized'
%                       and 'Pearson'.
%
%   For more details on the various residual types, please see the help for
%   the residuals method.
%  
%   H = plotResiduals(...) returns a handle to the lines or patches in the
%   plot.
%  
%   The PLOTTYPE or RTYPE arguments can be followed by parameter/value
%   pairs to specify additional properties of the primary line in the plot.
%   For example, plotResiduals(M,'fitted','Marker','s') uses a square
%   marker.
%     
%   For many of these plots, the data cursor tool in the figure window will
%   display the X and Y values for any data point, along with the
%   observation name or number.
%
%   Example: Fit a model and make a probability plot of the residuals.
%            There are two outliers in the upper right of the plot, which
%            you can examine using the data cursor in the plot figure.
%      load carsmall
%      obsname = strcat(num2str((1:100)'),{' '},Model);
%      ds = dataset(MPG,Weight,Model_Year,'ObsNames',obsname);
%      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
%      plotResiduals(lme,'probability')
%
%   See also RESIDUALS, FITTED.

        % (1) Set default values for positional parameters.
        if nargin < 2
            plottype = 'histogram';
        end
        
        % (2) Set parameter defaults.
        dfltResidualType = 'Raw';                

        % (3) Optional parameter names and their default values.
        paramNames =   {'ResidualType'};
        paramDflts = {dfltResidualType};
           
        % (4) Parse optional parameter name/value pairs.
        [residualtype,~,args] ...
            = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});  
        
        % (5) args may contain plot options. Save them in varargin.
        varargin = args;
        
        % (6) Validate residualtype.
        residualtype = model.validateResidualType(residualtype);
        
        % (7) Ensure that plot options given in varargin are sensible.
        internal.stats.plotargchk(varargin{:});

        % (8) Call helper method.
        f = classreg.regr.modelutils.plotResiduals(model,plottype,'ResidualType',residualtype,varargin{:});
        if nargout > 0
            hout = f;
        end

    end % end of plotResiduals.

    function stats = anova(model,varargin)
%ANOVA Perform hypothesis tests on fixed effect terms.
%   STATS = ANOVA(LME) tests the significance of each fixed effect term in
%   the linear mixed effects model LME and returns a table array STATS.
%   Each fixed effect term reported in STATS is either a continuous
%   variable, a grouping variable or an interaction between two or more
%   variables (continuous or grouping). For each fixed effect term, ANOVA
%   performs an F test (marginal test) that all coefficients representing
%   the fixed effect term are zero. There is one row in STATS for each
%   fixed effect term and the following columns:
%
%       Term        Name of the fixed effect term
%       FStat       F-statistic for the term
%       DF1         Numerator degrees of freedom for the F-statistic
%       DF2         Denominator degrees of freedom for the F-statistic
%       pValue      p-value for the term
%
%   To obtain tests for the Type III hypotheses, set the 'DummyVarCoding' 
%   to 'effects' when fitting the model.
%
%   STATS = ANOVA(LME,PARAM1,VALUE1,...) accepts additional parameter
%   name/value pairs:
%  
%           'Name'     'Value'
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate denominator degrees of freedom DF2 for
%                      the F-statistics reported in STATS. Options are
%                      'Satterthwaite','Residual' and 'None'. If 'DFMethod'
%                      is 'Satterthwaite', a Satterthwaite approximation is
%                      used to compute DF2. If 'DFMethod' is 'Residual',
%                      the DF2 values are assumed to be constant and equal
%                      to (N-P) where N is the number of observations and P
%                      is the number of fixed effects. If 'DFMethod' is
%                      'None', then all DF2 values are set to infinity.
%                      Default is 'Residual'.
%
%   Example: Fit a model with two fixed-effect predictors and compute an
%            anova table. Note that the p-values for the intercept and
%            Weight terms are the same as in the coefficients table. The
%            p-value for the Cylinders term measures the combined
%            significance of both Cylinders coefficients.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year,Cylinders);
%      ds.Cylinders = nominal(ds.Cylinders);
%      lme = fitlme(ds,'MPG ~ Weight + Cylinders + (1|Model_Year)')
%      anova(lme)
%
%   See also coefTest, coefCI, fixedEffects.

        % TODO: Add support for other types of sum of squares. 
        % (1) Default parameter values.
        dfltDFMethod = 'Residual';

        % (2) Optional parameter names and their default values.
        paramNames =   {'DFMethod'};
        paramDflts = {dfltDFMethod};
           
        % (3) Parse optional parameter name/value pairs.
        dfmethod ...
            = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % (4) Validate optional parameter values.
        dfmethod = model.validateDFMethod(dfmethod);

        % (5) What are the fixed effects term names?       
        termnames = model.Formula.FELinearFormula.TermNames;
        if size(termnames,1) == 1
            termnames = termnames';
        end
        nterms = length(model.Formula.FELinearFormula.TermNames);        
        
        % (6) Loop over the fixed effects terms. Create a contrast matrix
        % for an F-test on fixed effects coefficients related to each term.
        %
        % Suppose termkcols is a 1-by-p logical vector such that
        % termkcols(j) = 1 if column j is related to term k. Then we would
        % like to perform the test: L*beta = 0 where L = I(termkcols,:) and
        % I is the p-by-p identity matrix.
       
            % (1) Initialize I, P, F, DF1 and DF2.
            I = eye(model.slme.p);
            P = zeros(nterms,1);
            F = zeros(nterms,1);
            DF1 = zeros(nterms,1);
            DF2 = zeros(nterms,1);
            
            for k = 1:nterms
                % (2) Which columns of X are involved in term k?
                termkcols = (model.FixedInfo.XCols2Terms == k);

                % (3) Extract rows of I, with index in termkcols.
                L = I(termkcols,:);
                e = zeros(size(L,1),1);
                [P(k),F(k),DF1(k),DF2(k)] = ...
                    betaFTest(model.slme,L,e,dfmethod);
            end
        
        % (7) Assemble the output table.
        stats = table(termnames,F,DF1,DF2,P,'VariableNames',{'Term','FStat','DF1',...
            'DF2','pValue'});
%       Term        Name of the fixed effects term
%       FStat       F-statistic for the term
%       DF1         Numerator degrees of freedom for the F-statistic
%       DF2         Denominator degrees of freedom for the F-statistic
%       pValue      p-value for the term
            
        % (8) Add title to the ANOVA table.
        % <entry key="Title_anova">ANOVA type III or marginal tests: DFMethod = ''{0}''</entry>
        ttl = getString(message('stats:LinearMixedModel:Title_anova',dfmethod));
        stats = classreg.regr.lmeutils.titleddataset(stats,ttl);

    end % end of anova.
                   
end % end of methods(Access='public').

% Subclasses must define the following public methods.
methods (Abstract,Access='public')

    yfit = fitted(model,varargin);
    res = residuals(model,varargin);
    table = compare(model,altmodel,varargin);
    Y = response(model);
    
end

end % end of LinearLikeMixedModel.

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
