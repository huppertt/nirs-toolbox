classdef (Sealed) RepeatedMeasuresModel < matlab.mixin.CustomDisplay & classreg.learning.internal.DisallowVectorOps
%RepeatedMeasuresModel - Repeated measures model.
%   RM = FITRM(T,MODELSPEC) fits the model specified by MODELSPEC to
%   data in the table T, and returns the RepeatedMeasuresModel RM.
%   T is a table containing the values of the response variables and
%   the between-subject factors to be used as predictors in the model.
%   MODELSPEC specifies the response variable names and model terms
%   as a string such as 'y1-y5 ~ x1 + x2 + x3*x4'.
%
%   RepeatedMeasuresModel methods:
%       anova       - Analysis of variance for between-subject effects.
%       coeftest    - Hypothesis test on coefficients.
%       margmean    - Estimated marginal means.
%       epsilon     - Epsilon adjustment for repeated measures anova.
%       grpstats    - Descriptive statistics by group.
%       manova      - Multivariate analysis of variance.
%       mauchly     - Mauchly's test of sphericity.
%       multcompare - Multiple comparisons of marginal means.
%       plot        - Plot data with optional grouping.
%       plotprofile - Plot expected marginal means with optional grouping.
%       predict     - Compute predicted values.
%       random      - Generate new random response values.
%       ranova     - Repeated measures analysis of variance.
%
%   RepeatedMeasuresModel properties:
%       BetweenDesign      - Design for between-subject factors.
%       BetweenModel       - Model for between-subject factors.
%       BetweenFactorNames - Names of between-subject factors.
%       ResponseNames      - Names of responses.
%       WithinDesign       - Design for within-subject factors.
%       WithinModel        - Model for within-subject factors.
%       WithinFactorNames  - Names of within-subject factors.
%       Coefficients       - Table of estimated coefficients.
%       Covariance         - Table of estimated response covariances.
%       DFE                - Degrees of freedom for error.
% 
%   See also FITRM.

%   Copyright 2013-2014 The MathWorks, Inc.

    properties(Hidden,Constant)
        SeparateMeans = 'separatemeans';
        OrthogonalContrasts = 'orthogonalcontrasts';
        MeanResponse = 'meanresponse';
    end
    properties(GetAccess=public, SetAccess=private)
%BetweenDesign - Design for between-subjects factors.
%   The BetweenDesign property is a table containing the values of the
%   between-subject factors and the repeated measures.
%
%   See also BetweenModel, WithinDesign, WithinModel.
        BetweenDesign
        
%DFE - Degrees of freedom for error.
%   The DFE property is the degrees of freedom for error. DFE is the number
%   of observations minus the number of estimated coefficients in the
%   between-subjects model.
%
%   See also BetweenModel.
        DFE
    end
    properties(Access=public)
%WithinModel - Model for within-subjects factors.
%   The WithinModel property is a text representation of the model as a
%   function of the between-subject factors. This text may be a formula
%   involving the within-subject factors, or it may be any of the following
%   special strings:
%
%      'separatemeans'         Simple model that fits a separate mean to
%                              each repeated measure (default).
%      'orthogonalcontrasts'   Valid only when the WithinDesign has a
%                              single numeric factor T. Responses are
%                              the average, the slope of centered T, and in
%                              general all orthogonal contrasts for a
%                              polynomial up to T^(R-1), where R is the
%                              number of rows in the within-subject design.
%
%   The WithinModel property is not used in fitting the model to data.
%   Instead, it is the default within-subject model used with the MANOVA
%   and PREDICT methods.
%
%   See also BetweenModel, BetweenDesign, WithinDesign.
        WithinModel = RepeatedMeasuresModel.SeparateMeans
    end
    properties(Dependent)
        %WithinDesign - Design for within-subjects factors.
        %   The WithinDesign property is a table containing the values of the
        %   within-subject factors.
        %
        %   See also WithinModel, BetweenModel, BetweenDesign.
        WithinDesign
    end
    properties(Dependent, SetAccess=private)
%BetweenModel - Model for between-subjects factors.
%   The BetweenModel property is a text representation of the model as a
%   function of the between-subject factors.
%
%   This text is part of a formula representation of the model supplied to
%   FITRM. It is the text that appears to the right of the tilde.
%
%   See also BetweenDesign, WithinDesign, WithinModel.
        BetweenModel
        
%ResponseNames - Names of the responses.
%   The ResponseName property is a cell array of strings containing the
%   names of the variables used as responses.
%
%   See also BetweenFactorNames, WithinFactorNames.
        ResponseNames
        
%WithinFactorNames - Names of the within-subject factors.
%   The WithinFactorNames property is a cell array of strings containing
%   the names of the variables used as within-subject factors.
%
%   See also BetweenFactorNames, ResponseNames.
        WithinFactorNames
        
%BetweenFactorNames - Names of the between-subject factors.
%   The BetweenFactorNames property is a cell array of strings containing
%   the names of the variables used as between-subject factors.
%
%   See also WithinFactorNames, ResponseNames.
        BetweenFactorNames
        
%Coefficients - Estimated coefficients.
%   The Coefficients property is a table containing the values of the
%   estimated coefficients for fitting the repeated measures as a function
%   of the terms in the between-subjects model.
%
%   The coefficients for a categorical term are defined using "effects"
%   coding. There is one coefficient for each level except the first. The
%   implied coefficient for the first level is the sum of the other
%   coefficients for the term. Use the MARGMEAN method to display marginal
%   means for all levels.
%
%   Example: Fit a repeated-measures model and create a matrix of
%            coefficient estimates from the Coefficients table.
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within);
%      R.Coefficients           % as table for easy viewing
%      B = R.Coefficients{:,:}  % as matrix for further calculations
%
%   See also Covariance.
        Coefficients

%Covariance - Estimated covariance.
%   The Covariance property is a table containing the estimated covariance
%   of the repeated measures. This covariance is computed around the mean
%   given by the fitted model.
%
%   Example: Fit a repeated-measures model and create a matrix of
%            covariance estimates from the Coefficients table.
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within);
%      R.Covariance           % as table for easy viewing
%      B = R.Covariance{:,:}  % as matrix for further calculations
%
%   See also Coefficients.
        Covariance

%DesignMatrix - Design matrix.
%   The DesignMatrix property is an R-by-T matrix where R is the number of
%   rows in the BetweenDesign. Row J contains the values of all the
%   BetweenModel terms for subject J. Any rows in the BetweenDesign that
%   contain missing (NaN) values are omitted from the design matrix.
%
%   Example: Fit a repeated-measures model and examine the design matrix.
%            It has one column for each term in the BetweenModel.
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%      disp(R.DesignMatrix)
%
%   See also Coefficients.
        DesignMatrix
    end
    properties(Access=private)
        ResponseColumns
        Mauchly
        Epsilon
        TermAverages
        IsCat
        VariableRange
        WithinDesign_
        Terms
        Missing
        X
        Y
        B
        Cov
        CoefNames
        CoefTerms
        TermNames
        Formula
    end
    
    
    % --- Public methods ---
    methods
        function tbl = grpstats(this,grp,stats)
%GRPSTATS Descriptive statistics by group
%   TBL = GRPSTATS(RM,G) computes the mean and variance for the data in the
%   repeated measures model RM, grouped by the factors G. G can be a factor
%   name or a cell array of factor names. The result TBL is a table giving
%   the mean and variance for each group.
%
%   TBL = GRPSTATS(RM,G,STATS) returns the statistics specified by STATS.
%   STATS can be a single function handle or name, or a cell array
%   containing multiple function handles or names. Names in STATS can be
%   chosen from among the following:
% 
%       'mean'     mean
%       'sem'      standard error of the mean
%       'numel'    count, or number of elements
%       'gname'    group name
%       'std'      standard deviation
%       'var'      variance
%       'min'      minimum
%       'max'      maximum
%       'range'    maximum - minimum
%       'meanci'   95% confidence interval for the mean
%       'predci'   95% prediction interval for a new observation
% 
%    Each function included in STATS must accept a vector of response
%    values for a single group, and compute descriptive statistics for it.
%    A function should typically return a value that has one row. A
%    function must return the same size output each time GRPSTATS calls it,
%    even if the input for some groups is empty.
%
%    GRPSTATS computes results separately for each group. The results do
%    not depend on the fitted repeated measures model. They are computed on
%    all available data, not omitting entire rows that contain NaN.
%
%    Example: Compute statistics by two between-subjects factors.
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%      T = grpstats(R,{'Group','Gender'})
%
%    See also PLOT.

            if nargin<3
                stats = {'mean' 'std'};
            elseif ~iscell(stats)
                stats = {stats};
            end
            grouped = ~(nargin<2 || isempty(grp));
            if grouped
                [grp,iswithin] = getgroup(this,grp);
            else
                grp = [];
            end
            
            % Stack the data to combine all responses into a y_ variable
            % and to introduce a time variable that represents the response
            % number.
            yname = genvarname('y_',this.BetweenDesign.Properties.VariableNames);
            timename = genvarname('time',this.BetweenDesign.Properties.VariableNames);
            tbl = this.BetweenDesign;
            dstacked = stack(tbl,this.ResponseColumns,'newdatavar',yname,'indexv',timename);
            
            if grouped
                % copy any required within-subjects factors
                ny = sum(this.ResponseColumns);
                nsubjects = size(tbl,1);
                wrows = repmat((1:ny)',nsubjects,1);
                dstacked = [dstacked, this.WithinDesign(wrows,grp(iswithin))];
            end
            
            % Run the regular grpstats function on this table
            tbl = grpstats(dstacked,grp,stats,'datavar',yname);
            
            % Fix the names to remove the y variable suffix and replace by
            % the statistics name. Note that the gname column is placed
            % earlier in the table, not among the other statistics.
            nstats = length(stats);
            statNames = cell(1,nstats);
            for j=1:nstats
                statNames{j} = char(stats{j});
            end
            statNames(strcmp('gname',statNames)) = [];
            nstats = length(statNames);
            vnames = tbl.Properties.VariableNames;
            vnames{end-nstats} = 'GroupCount';
            vnames(end-nstats+1:end) = statNames;
            tbl.Properties.VariableNames = vnames;
            tbl.Properties.RowNames = {};
        end

        function [tbl,A,C,D] = manova(this,varargin)
%MANOVA Multivariate analysis of variance.
%   TBL = MANOVA(RM) produces a table TBL containing the multivariate
%   analysis of variance (manova) for the repeated measures model RM. This
%   manova includes all terms in the between-subjects model. The
%   multivariate response for each observation (subject) is the vector of
%   repeated measures.
%
%   TBL = MANOVA(RM,NAME1,VAL1,...) also specifies one or more of the
%   following name/value pairs:
% 
%        'WithinModel'  A model MODEL specifying the within-subjects
%                       hypothesis test. Default is the model specified
%                       when creating RM. MODEL may be any of the following:
%
%            'separatemeans'   Compute a separate mean for each group, and
%                              test for equality among the means.
%            C                 An R-by-NC matrix specifying NC contrasts
%                              among the R repeated measures.
%            MODELSPEC         A string that defines a model specification
%                              in the within-subject factors. See FITRM for
%                              a description of a model specification. The
%                              terms in that model define the columns of
%                              the contrast matrix C.
%
%        'By'           A single between-subjects factor B. The output TBL
%                       contains a separate hypothesis test on the
%                       within-subjects model for each value of B.
%
%   [TBL,A,C,D] = MANOVA(...) also returns arrays A, C, and D. The
%   hypotheses tested by this function have the form A*B*C=D, where B is
%   the coefficient matrix, A is determined by the between-subjects model,
%   C is determined by the within-subjects model, and D is zero. The
%   estimated value of B is given by RM.Coefficients. The outputs are the
%   matrices used in this test. If TBL contains multiple hypothesis tests,
%   the A and C outputs may be cell arrays of matrices.
%
%   Example:
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%      T = manova(R)
%
%   See also ANOVA, RANOVA, FITRM, COEFTEST.
            okargs =   {'By' 'WithinModel'};
            defaults = {''   this.WithinModel};
            [by,withinmodel] = internal.stats.parseArgs(okargs,defaults,varargin{:});
            if strcmp(withinmodel,RepeatedMeasuresModel.OrthogonalContrasts)
                error(message('stats:fitrm:NoOrthogonalManova'));
            end
            
            % Get C matrix and within-subject names
            [C,wnames] = makeTestC(this,withinmodel,true,this.WithinDesign,true);
            
            % Get A matrix and between-subject names
            if isempty(by)
                % separate tests for each between-subjects term
                [A,bnames] = makeTestA(this);
            else
                % separate tests for each level of the BY variable
                [A,bnames] = makeTestABy(this,by);
            end
            
            % Get manova table
            Beta = this.Coefficients{:,:};
            SSE = this.DFE * this.Cov;
            tbl = manovastats(this.X,A,Beta,C,0,SSE,bnames,wnames);
            tbl.Properties.Description = getString(message('stats:fitrm:TableDescrManova'));
            D = 0;
        end
        
        function tbl = coeftest(this,A,C,D)
%COEFTEST Linear hypothesis test on coefficients.
%   TBL = COEFTEST(RM,A,C,D) produces a table TBL containing the
%   multivariate analysis of variance (manova) for the repeated measures
%   model RM. This test is defined as
%
%         A*B*C = D
%
%   where B is the matrix of coefficients in the repeated measures model. A
%   and C are numeric matrices of the proper size for this multiplication.
%   D is a scalar or numeric matrix of the proper size. Default is D=0.
%
%   Example: Test that the coefficients of all terms in the between-
%            subjects model are the same for the first and last repeated
%            measurement variable.
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%      coeftest(R,eye(8),[1 0 0 0 0 0 0 -1]')
%
%   See also ANOVA, RANOVA, FITRM.
            nx = size(this.X,2);
            ny = size(this.Y,2);
            checkMatrix('A',A,[],nx);
            checkMatrix('C',C,ny,[]);
            if nargin<4
                D = 0;
            else
                checkMatrix('D',D,size(A,1),size(C,2));
            end

            Beta = this.Coefficients{:,:};
            SSE = this.DFE * this.Cov;
            tbl = fourstats(this.X,A,Beta,C,D,SSE);
        end
        
        function tbl = mauchly(this,C)
%MAUCHLY Mauchly's test of sphericity.
%   TBL = MAUCHLY(RM) performs Mauchly's test of sphericity for the
%   repeated measures model RM. The output TBL is a single-row table
%   containing four variables:
%
%       W       Mauchly's W statistic
%       ChiStat Chi-squared test statistic
%       DF      Degrees of freedom
%       pValue  p-value
%
%   TBL = MAUCHLY(RM,C) performs Mauchly's test based on the contrast
%   matrix C. The default value of C is the Q factor in a QR decomposition
%   of M, where M is defined so that Y*M is the difference between all
%   successive pairs of columns of the repeated measures matrix Y.
%
%   For a repeated measures model with responses Y1, Y2, Y3, ...,
%   sphericity means that all pairwise differences Y1-Y2, Y1-Y3, Y2-Y3, ...
%   have the same theoretical variance.
%
%   The regular repeated measures anova p-values are valid under the
%   assumption of compound symmetry, which is that all Y responses have the
%   same variance and all pairs have a common correlation. Compound
%   symmetry implies sphericity. If the assumption of compound symmetry is
%   not true, then the degrees of freedom for the repeated measures anova
%   test must be adjusted by a factor known as epsilon, and the p-value
%   must be computed using the adjusted values.
%
%   The EPSILON method returns a table of epsilon adjustment values. The
%   table produced by the RANOVA method contains p-values based on each
%   epsilon value.
%
%    Example:
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%      T = mauchly(R)
%
%   See also RANOVA, EPSILON, FITRM.
            if nargin<2
                tbl = this.Mauchly;
            else
                if ~iscell(C)
                    C = {C};
                end
                Xmat = this.X;
                ny = size(this.Y,2);
                for j=1:numel(C)
                    Cj = C{j};
                    checkMatrix('C',Cj,ny,[]);
                    tbl(j,:) = mauchlyTest(this.Cov,size(Xmat,1),rank(Xmat),Cj);
                end
            end
        end
        
        function tbl = epsilon(this,C)
%EPSILON Epsilon adjustment for repeated measures anova.
%   TBL = EPSILON(RM) produces a table TBL of epsilon adjustment factors
%   for the repeated measures model RM. The table has a single variable
%   named epsilon, and four rows:
%
%      Uncorrected            No correction, epsilon=1
%      Greenhouse-Geisser     Greenhouse-Geisser approximation
%      Huynh-Feldt            Huynh-Feldt approximation
%      Lower bound            Lower bound on the true p-value
%
%   TBL = EPSILON(RM,C) produces a table of adjustment factors for the test
%   based on the contrast matrix C. The default value of C the Q factor in
%   a QR decomposition of M, where M is defined so that Y*M is the
%   difference between all successive pairs of columns of the repeated
%   measures matrix Y.
%
%   For a repeated measures model with responses Y1, Y2, Y3, ...,
%   sphericity means that all pairwise differences Y1-Y2, Y1-Y3, Y2-Y3, ...
%   have the same theoretical variance.
%
%   The regular repeated measures anova p-values are valid under the
%   assumption of compound symmetry, which is that all Y responses have the
%   same variance and all pairs have a common correlation. Compound
%   symmetry implies sphericity. If the assumption of compound symmetry is
%   not true, then the degrees of freedom for the repeated measures anova
%   test must be adjusted by a factor known as epsilon, and the p-value
%   must be computed using the adjusted values.
%
%   The MAUCHLY method tests for sphericity. The RANOVA method contains
%   p-values based on each epsilon value.
%
%    Example:
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%      T = epsilon(R)
%
%   See also RANOVA, EPSILON, FITRM.
            if nargin<2
                tbl = this.Epsilon;
            else
                if ~iscell(C)
                    C = {C};
                end
                Xmat = this.X;
                ny = size(this.Y,2);
                for j=1:numel(C)
                   Cj = C{j};
                   checkMatrix('C',Cj,ny,[]);
                   [~,tbl(j,:)] = mauchlyTest(this.Cov,size(Xmat,1),rank(Xmat),Cj);
                end
            end
        end
        
        function [tbl,Q] = anova(this,varargin)
%ANOVA Analysis of variance for between-subject effects.
%   TBL = ANOVA(RM) produces a table TBL containing the analysis of
%   variance (anova) for the repeated measures model RM. This anova
%   includes all terms in the between-subjects model. The response is
%   calculated by averaging across values of the within-subjects factors.
%
%   TBL = ANOVA(RM,'WithinModel',WM) performs anova using the response or
%   responses specified by the within-subject model WM. WM may be any of
%   the following:
%
%      'meanresponse'          Response is averaged across within-subject
%                              design (this is the default).
%      'orthogonalcontrasts'   Valid only when the within-subject design
%                              has a single numeric factor T. Responses are
%                              the average, the slope of centered T, and in
%                              general all orthogonal contrasts for a
%                              polynomial up to T^(R-1), where R is the
%                              number of rows in the within-subject design.
%      C                       An R-by-NC matrix specifying NC contrasts
%                              among the R repeated measures. If Y
%                              represents a matrix of repeated measures,
%                              the output TBL contains a separate anova for
%                              each column of Y*C.
%      MODELSPEC               A string that defines a model specification
%                              in the within-subject factors. See FITRM for
%                              a description of a model specification.
%                              The terms in that model define the columns
%                              of the contrast matrix C.
%
%   This ANOVA table contains blocks of rows with separate univariate ANOVA
%   results for each response. Orthogonal contrasts for T are computed
%   using the Q factor of a QR factorization of the Vandermonde matrix.
%
%    Example:
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%
%      % anova on average across repeated measures
%      T1 = anova(R)
%
%      % separate anova for two successive time differences
%      C = [1 -1 0 0 0 0 0 0;0 1 -1 0 0 0 0 0]'
%      T2 = anova(R,'WithinModel',C)
%
%   See also RANOVA, MANOVA, FITRM, QR, VANDER.

% The following is likely to change in a future release:
%   [TBL,C] = ANOVA(...) returns the contrast matrix C that is used to
%   generate the response values for each anova.

            okargs =   {'WithinModel'};
            defaults = {RepeatedMeasuresModel.MeanResponse};
            wm = internal.stats.parseArgs(okargs,defaults,varargin{:});
            if isempty(wm)
                wm = RepeatedMeasuresModel.MeanResponse;
            end

            ny = size(this.Y,2);
            d = this.BetweenDesign(~this.Missing,:);
            yname = genvarname('y',this.BetweenDesign.Properties.VariableNames);
            newformula = sprintf('%s ~ %s',yname,this.Formula.LinearPredictor);
            if isnumeric(wm)
                % We have a matrix to be used to define contrasts among the
                % repeated measures, one contrast per column
                checkMatrix('WithinModel',wm,ny,[]);
                Q = wm;
                contrastNames = textscan(sprintf('Contrast%d\n',1:size(Q,2)),'%s');
                contrastNames = contrastNames{1};
            else
                % We have a definition of a model, and the terms in that
                % model are to be used to define contrasts
                [C,contrastNames] = makeTestC(this,wm,false);
                Q = bsxfun(@rdivide,C,sqrt(sum(C.^2,1)));
            end
            tbl = contrastanova(d,newformula,this.Y*Q,contrastNames,false,yname);
            tbl.Properties.Description = getString(message('stats:fitrm:TableDescrAnova'));
        end

        function [tbl,A,C,D] = ranova(this,varargin)
%RANOVA Repeated measures analysis of variance.
%   TBL = RANOVA(RM) produces a table TBL containing the repeated measures
%   analysis of variance (anova) for the repeated measures model RM. This
%   anova includes all terms in the between-subjects model, each
%   interacting with a single Time term that represents all differences
%   across the within-subject factors.
%
%   TBL = RANOVA(RM,'WithinModel',WM) performs anova using the response or
%   responses specified by the within-subject model WM. WM may be any of
%   the following:
%
%      'separatemeans'         Test differences between separate means
%                              (this is the default).
%      C                       An R-by-NC matrix specifying NC contrasts
%                              among the R repeated measures. If Y
%                              represents a matrix of repeated measures,
%                              the test is that the means of Y*C are zero.
%      MODELSPEC               A string that defines a model specification
%                              in the within-subject factors. See FITRM for
%                              a description of a model specification.
%                              The terms in that model define the columns
%                              of the contrast matrix C.
%
%   [TBL,A,C,D] = RANOVA(...) also returns arrays A, C, and D. The
%   hypotheses tested by this function have the form A*B*C=D, where B is
%   the coefficient matrix, A is determined by the between-subjects model,
%   C is determined by the within-subjects model, and D is zero. The
%   estimated value of B is given by RM.Coefficients. The outputs are the
%   matrices used in this test. If TBL contains multiple hypothesis tests,
%   the A and C outputs may be cell arrays of matrices.
%
%   The pValue in TBL is a measure of the significance of each term. This
%   p-value is accurate for data with "compound symmetry." In case the
%   compound symmetry assumption does not hold, the output TBL contains
%   three other variables with p-values intended to adjust for failure to
%   satisfy the assumption:
%
%      pValueGG   Greenhouse-Geisser adjustment
%      pValueHF   Huynh-Feldt adjustment
%      pValueLB   Lower bound adjustment
%
%   For a repeated measures model with responses Y1, Y2, Y3, ...,
%   sphericity means that all pairwise differences Y1-Y2, Y1-Y3, Y2-Y3, ...
%   have the same theoretical variance.
%
%   The regular repeated measures anova p-values are valid under the
%   assumption of compound symmetry, which is that all Y responses have the
%   same variance and all pairs have a common correlation. Compound
%   symmetry implies sphericity. If the assumption of compound symmetry is
%   not true, then the degrees of freedom for the repeated measures anova
%   test must be adjusted by a factor known as epsilon, and the p-value
%   must be computed using the adjusted values.
%
%   The MAUCHLY method tests for sphericity. The EPSILON method returns a
%   table of epsilon adjustment values.
%
%    Example:
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%      T = ranova(R)
%
%   See also MAUCHLY, EPSILON, ANOVA, MANOVA, FITRM.
            okargs =   {'WithinModel'};
            defaults = {RepeatedMeasuresModel.SeparateMeans};
            wm = internal.stats.parseArgs(okargs,defaults,varargin{:});
            if isempty(wm)
                wm = RepeatedMeasuresModel.SeparateMeans;
            elseif isnumeric(wm)
                checkMatrix('WithinModel',wm,size(this.Y,2),[]);
            end

            if isempty(this.WithinDesign) || size(this.WithinDesign,2)>1
                timename = genvarname('Time',this.BetweenDesign.Properties.VariableNames);
            else
                timename = this.WithinDesign.Properties.VariableNames{1};
            end

            % Get A and C matrices for test of the form   A*B*C=D
            [A,Anames] = makeTestA(this);
            [C,Cnames] = makeTestC(this,wm,true,[],true);
            D = 0;
            
            if iscell(C)
                dscell = cell(size(C));
                nc = numel(C);
                na = numel(Anames);
                rownames = cell(nc*(na+1),1);
                baserow = 0;
                for j=1:numel(C)
                    tblj = ranovastats(this,C{j},A,Anames,D,Cnames{j});
                    rownames(baserow+(1:na+1),:) = tblj.Properties.RowNames;
                    dscell{j} = tblj;
                    baserow = baserow+na+1;
                end
                tbl = vertcat(dscell{:});
                tbl.Properties.RowNames = rownames;
            else
                tbl = ranovastats(this,C,A,Anames,D,timename);
            end
        end

        function hout = plot(this,varargin)
%PLOT Plot data with optional grouping.
%     PLOT(RM) plots the measurements in the RepeatedMeasuresModel RM
%     for each subject as a function of time. If there is a single
%     numeric within-subjects factor, its values are used as the time
%     values. Otherwise the time value is taken as 1:NY where NY is
%     the number of repeated measurements.
% 
%     PLOT(RM,NAME1,VAL1,...) also specifies one or more of the
%     following name/value pairs:
% 
%        'Group'     The name of a between-subjects factor, or a cell
%                    array of multiple names. The plotted lines are
%                    grouped according to the factor values.
%        'Marker'    A cell array of strings specifying the marker
%                    to be used for each group.
%        'Color'     The color for each group, specified as a cell
%                    array of strings or as the rows of a three-
%                    column RGB matrix.
%        'LineStyle' A cell array of strings specifying the line
%                    style for each group.
%
%     H = PLOT(RM,...) also returns handles H to the plotted lines.
%
%    Example: Plot with Group coded by color and Gender by line type.
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%      plot(R,'group',{'Group' 'Gender'},'Color','rrbbgg',...
%             'LineStyle',{'-' ':' '-' ':' '-' ':'},'Marker','.')
%
%    See also PLOTPROFILE.
            markers = {'s','o','*','x','+','d','^','v','>','<','p','h'};
            okargs =   {'Group' 'Marker' 'Color', 'LineStyle'};
            defaults = {''      markers  ''       {'-'}};
            [group,markers,cmap,styles] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});
            w = this.WithinDesign;
            if ~isempty(w) && size(w,2)==1 ...
                           && varfun(@(x)isnumeric(x), w,'OutputFormat','uniform')
                x = w{:,1};
                xticks = [];
            else
                x = (1:size(this.Y,2));
                xticks = x;
            end
            
            grouped = ~isempty(group);
            [ngroups,grpidx,grpname] = makegroup(group,this);
            [cmap,markers,styles] = regularizePlotArgs(cmap,markers,styles,ngroups);
            newplot;
            h = [];
            hleg = [];
            for j=1:ngroups
                idx = grpidx==j;
                jcolor = 1 + mod(j-1,size(cmap,1));
                jmarker = 1 + mod(j-1,numel(markers));
                jstyle = 1 + mod(j-1,numel(styles));
                hj = line(x,this.Y(idx,:)','Color',cmap(jcolor,:),...
                    'Marker',markers{jmarker}, 'LineStyle',styles{jstyle});
                h = [h; hj];          %#ok<AGROW>
                hleg = [hleg; hj(1)]; %#ok<AGROW>
            end
            if grouped
                legend(hleg,grpname,'Location','best')
            end
            xlim = get(gca,'XLim');
            dx = diff(xlim)/20;
            set(gca,'XLim',[xlim(1)-dx,xlim(2)+dx])
            if ~isempty(xticks)
                set(gca,'XTick',xticks);
            end
            if nargout>0
                hout = h;
            end
        end

        function hout = plotprofile(this,x,varargin)
%PLOTPROFILE Plot expected marginal means with optional grouping.
%     PLOTPROFILE(RM,X) plots the expected marginal means computed from the
%     RepeatedMeasuresModel RM as a function of the variable X. X is a text
%     string that gives the name of a between- or within-subjects factor.
% 
%     PLOTPROFILE(RM,NAME1,VAL1,...) also specifies one or more of the
%     following name/value pairs:
% 
%        'Group'     The name of a factor, or a cell array of multiple
%                    names. The plotted lines are grouped according to
%                    the factor values.
%        'Marker'    A cell array of strings specifying the marker
%                    to be used for each group.
%        'Color'     The color for each group, specified as a cell
%                    array of strings or as the rows of a three-
%                    column RGB matrix.
%        'LineStyle' A cell array of strings specifying the line
%                    style for each group.
% 
%     H = PLOTPROFILE(RM,...) also returns handles H to the plotted lines.
%
%    Example:
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%      ax1 = subplot(1,2,1);
%      plotprofile(R,'Group')
%      ax2 = subplot(1,2,2);
%      plotprofile(R,'Gender')
%      linkaxes([ax1 ax2],'y')
%
%    See also PLOT.
            okargs =   {'Group' 'Marker' 'Color', 'LineStyle'};
            defaults = {{}      {'o'}    ''       {'-'}};
            [group,markers,cmap,styles] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});

            if ~internal.stats.isString(x,false) || ~isVariable(this,x)
                error(message('stats:fitrm:ValueMustBeFactor','X'));
            end
            if iscell(group)
                group = group(:)';
                vars = [{x} group];
            else
                vars = {x group};
            end
            
            % Compute means sorted first by group
            m = margmean(this,vars);
            m = sortrows(m,fliplr(vars));

            % Plot against x data if numeric, or category number otherwise
            xvals = m.(x);
            if ischar(xvals)
                xvals = cellstr(xvals);
            end
            [ux,~,uidx] = unique(xvals,'stable');
            nx = length(ux);
            xloc = 1:nx;
            if ~isnumeric(xvals)
                xvals = uidx;
            end
            
            if isempty(group)
                % Plot without grouping
                h = plot(xvals,m.Mean,'o-');
            else
                % Separate lines for each group
                nrows = size(m,1);
                R = nx;
                C = nrows/R;
                means = reshape(m.Mean,[R,C]);
                h = plot(xvals(1:nx),means,'o-');
                if iscell(group)
                    gvals = cell(size(group));
                    for j=1:numel(group)
                        gvals{j} = m.(group{j});
                    end
                    [~,~,gcell] = internal.stats.mgrp2idx(fliplr(gvals));
                    gcell = fliplr(gcell);
                    GN = gcell(:,1);
                    for j=2:size(gcell,2);
                        GN = strcat(GN,',',gcell(:,j));
                    end
                    varnames = sprintf('%s,',group{:});
                    varnames(end) = [];
                else
                    gvals = m.(group);
                    [~,GN] = grp2idx(gvals);
                    varnames = group;
                end
                for j=1:length(h)
                    set(h(j),'DisplayName',sprintf('%s=%s',varnames,GN{j}))
                end
                legend('Location','best');
            end
            [cmap,markers,styles] = regularizePlotArgs(cmap,markers,styles,length(h));
            for j=1:length(h)
                jcolor = 1 + mod(j-1,size(cmap,1));
                jmarker = 1 + mod(j-1,numel(markers));
                jstyle = 1 + mod(j-1,numel(styles));
                set(h(j),'Color',cmap(jcolor,:),...
                    'Marker',markers{jmarker}, 'LineStyle',styles{jstyle});
            end
           
            if isa(ux,'categorical')
                set(gca,'XTick',xloc,'XTickLabel',char(ux));
            end
            
            % Usually the x values are small integers, so make sure the
            % default limits push them into the interior
            xmin = min(xvals);
            xmax = max(xvals);
            dx = (xmax-xmin)/20;
            set(gca,'XLim',[xmin-dx, xmax+dx])
            xlabel(x);
            ylabel(getString(message('stats:fitrm:EstimatedMarginalMeans')));

            if nargout>0
                hout = h;
            end
        end

        function ds = margmean(this,grp,varargin)
%MARGMEAN Estimated marginal means.
%   TBL = MARGMEAN(RM,VARS) produces a table TBL containing the estimated
%   marginal means for variables VARS in the repeated measures model RM.
%   VARS is the name of a between- or within-subjects factor in RM, or a
%   cell array of multiple names. Each between-subjects factor must be
%   categorical.
%
%   The output table TBL has one row for each combination of the variables
%   specified as VARS. It has one column for each of those variables, a
%   Mean column containing the estimated marginal means, a StdErr column
%   containing their standard errors, and Lower and Upper columns define
%   95% confidence intervals for the true (population) means.
%
%   TBL = MARGMEAN(RM,VARS,'Alpha',ALPHA) accepts a scalar 0<ALPHA<1 to
%   specify 100*(1-ALPHA)% confidence intervals. Default is ALPHA=0.05 for
%   95% confidence.
%
%    Example:
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ')
%      T = margmean(R,{'Group' 'Time'})
%      T.Properties.Description
%
%   See also PLOTPROFILE, FITRM, MULTCOMPARE.

            okargs =   {'Alpha'};
            defaults = {0.05};
            alpha = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});
            checkAlpha(alpha);

            % Find grouping variables. Separate the Between and Within
            % variables, and consider only Between for now.
            [grp,iswithin] = getgroup(this,grp,true);

            % Get full factorial and design matrix for the Between factors
            [A,Bff,~,cap] = emmMakeA(this,grp,iswithin);
            
            % Get full factorial and design matrix for the Within factors
            [C,Wff] = emmMakeC(this,grp,iswithin);
            
            % Compute the emm and its standard errors
            [mn,se] = meanCalculate(this,A,C);
            
            % Package up as a table
            ds = emmMakeDataset(Bff,Wff,grp,mn(:),se(:),this.DFE,alpha);
            ds.Properties.Description = sprintf('%s\n%s',...
                getString(message('stats:fitrm:EstimatedMarginalMeans')),cap);
        end

        function [ypred,yci] = predict(this,ds,varargin)
%PREDICT Compute predicted values.
%    YPRED = PREDICT(RM,T) computes a matrix YPRED of predicted values from
%    the repeated measures model RM using predictor variables taken from
%    the table T. T must contain all of the between-subject factors
%    used to create RM. The output YPRED is an N-by-R matrix, where N is
%    the number of rows in T and R is the number of repeated measures
%    in RM. The default for T is the table used to create RM.
%
%    YPRED = predict(RM,T,PARAM1,VAL1,...) specifies one or more
%    of the following name/value pairs: 
%
%       'Alpha'        Specifies the confidence level as 100*(1-ALPHA)%
%                      (default 0.05 for 95% confidence).
%       'WithinModel'  A model WMODEL specifying the model for the within-
%                      subjects factors. Default is the model specified
%                      when creating RM. WMODEL may be any of the following:
%            'separatemeans'
%                    Compute a separate mean for each group
%            'orthogonalcontrasts'
%                    Valid only when then the within-subject design
%                    consists of a single numeric factor. If that factor
%                    is Time, specifies a model consisting of orthogonal
%                    polynomials up to order Time.^(R-1).
%            WMODELSPEC
%                    Test string specifying the model terms for the within-
%                    subjects model.
%       'WithinDesign' A design for within-subject factors. It provides the
%                      values of the within-subject factors in the same
%                      form as the RM.WithinDesign property -- either as a
%                      vector, matrix, or table.
%
%    [YPRED,YCI] = PREDICT(...) also returns the N-by-R-by-2 array YCI
%    containing 95% confidence intervals for the predicted values. These are
%    non-simultaneous intervals for predicting the mean response at the
%    specified predictor values. For predicted value YPRED(I,J), the
%    lower limit of the interval is YCI(I,J,1) and the upper limit is
%    YCI(I,J,2).
%
%    Example:
%      % Predict at intermediate times
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group')
%      t = linspace(0,8)';
%      Y = predict(R,between([1 11 21],:), ...
%                  'WithinModel','Time^2','WithinDesign',t);
%
%      % Overlay onto profile plot
%      plotprofile(R,'Time','Group','Group')
%      hold on; plot(t,Y,'LineStyle',':'); hold off
%
%    See also FITRM, RANDOM.
            okargs =   {'Alpha' 'WithinDesign'     'WithinModel'};
            defaults = {0.05    this.WithinDesign  this.WithinModel};
            [alpha,withindesign,withinmodel,setflag] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});
            checkAlpha(alpha);
            if isa(withindesign,'dataset')
                withindesign = dataset2table(withindesign);
            elseif isrow(withindesign) && ~istable(withindesign) && numel(this.ResponseNames)>1
                withindesign = withindesign(:);
            end
            
            % Compute design matrix for input data and calculated predicted
            % values for all responses
            if nargin<2
                A = this.X;
            else
                if isa(ds,'dataset')
                    ds = dataset2table(ds);
                elseif ~isa(ds,'table')
                    error(message('stats:fitrm:TableRequired2'));
                end
                [terms,iscat,vrange] = dataset2terms(this,ds);
                A = classreg.regr.modelutils.designmatrix(ds,...
                    'Model',terms, ...
                    'DummyVarCoding','effects', ...
                    'CategoricalVars',iscat, ...
                    'CategoricalLevels',vrange);
            end

            if isempty(withindesign) || strcmp(withinmodel,RepeatedMeasuresModel.SeparateMeans)
                if ~isempty(withindesign) && setflag.WithinDesign
                    % Design is specified but model does not use it
                    warning(message('stats:fitrm:WithinDesignIgnored'));
                end
                C = eye(size(this.B,2));
            else
                % If a within-subjects design is given, predict using that.
                % Check the design and get the old and new time vectors.
                [withindesign,time] = verifyWithinDesign(this,withinmodel,withindesign);
                switch(withinmodel)
                    case RepeatedMeasuresModel.OrthogonalContrasts
                        % Compute the terms for the original design
                        W = fliplr(vander(time(:)));
                        [Q,R] = qr(W,0);

                        % Compute the terms for the new design points
                        ny = size(W,2);
                        V = ones(numel(withindesign),ny);
                        for j=2:ny
                            V(:,j) = withindesign.^(j-1);
                        end
                    otherwise
                        W = makeTestC(this,withinmodel,false,this.WithinDesign);
                        [Q,R] = qr(W,0);
                        V = makeTestC(this,withinmodel,false,withindesign);
                end
                
                % Compute the coefficients from the original design and
                % apply them to the new points
                C = Q/R'*V';
            end
            [ypred,yse] = meanCalculate(this,A,C);
            ypred = ypred';
            yse = yse';
            width = yse * -tinv(alpha/2,this.DFE);
            yci = cat(3,ypred-width,ypred+width);
        end

        function ysim = random(this,varargin)
%RANDOM Generate new random response values.
%    YNEW = RANDOM(RM,T) generates a matrix YNEW of random response
%    values from the repeated measures model RM using predictor variables
%    taken from the table T. T must contain all of the between-subject
%    factors used to create RM. The output YSIM is an N-by-R matrix, where
%    N is the number of rows in T and R is the number of repeated measures
%    in RM. YSIM is computed by creating predicted values and adding
%    new random noise values. For each row, the noise has a multivariate
%    normal distribution with covariance RM.Covariance.
%
%    The default for T is the table used to create RM.
%
%    Example:
%      % Predict at intermediate times
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%      Y = random(R);
%
%    See also FITRM, PREDICT.
            if nargin>2
                error(message('MATLAB:TooManyInputs'))
            end
            ymean = predict(this,varargin{:});
            ysim = mvnrnd(ymean,this.Cov);
        end
        
        function tbl = multcompare(this,var,varargin)
%MULTCOMPARE Multiple comparison of estimated marginal means.
%   TBL = MULTCOMPARE(RM,VAR) performs multiple comparisons of the
%   estimated marginal means based on the variable VAR in the repeated
%   measures model RM. VAR must be the name of a within-subjects factor or
%   a categorical between-subjects factor in RM. The output table TBL
%   provides a comparison between the estimated marginal means for each
%   distinct ordered pair of levels of VAR. The comparison includes the
%   estimated difference, its standard error, a p-value for a test that the
%   difference is zero, and lower and upper limits of simultaneous 95%
%   confidence intervals for the true difference.
%
%   TBL = MULTCOMPARE(RM,VAR, 'PARAM1',val1, 'PARAM2',val2,...) specifies
%   one or more of the following name/value pairs:
%   
%     'Alpha'       Specifies the confidence level as 100*(1-ALPHA)%
%                   (default 0.05).
%
%     'By'          A single factor B. The comparison between levels of VAR
%                   is done separately for each value of B.
%
%     'ComparisonType'  The type of critical value to use.  Choices are
%                   'tukey-kramer' (default), 'dunn-sidak', 'bonferroni',
%                   'scheffe', or 'lsd'. The 'lsd' option stands for "least
%                   significant difference" and uses plain t-tests; it
%                   provides no protection against the multiple comparison
%                   problem.
%
%    Example:
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%      T = multcompare(R,'Group')
%
%    See also MARGMEAN.
            okargs =   {'Alpha' 'ComparisonType' 'By'};
            defaults = {0.05    'tukey-kramer'   {}};
            [alpha,ctype,by] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});
            checkAlpha(alpha);
            
            % Find grouping variables.
            [var,iswithin] = getgroup(this,var,true);
            if ~isscalar(iswithin)
                error(message('stats:fitrm:ValueMustBeVar','VAR'));
            end
            isby = ~isempty(by);
            if isby
                [by,bywithin] = getgroup(this,by,true);
                if ~isscalar(bywithin)
                    error(message('stats:fitrm:ValueMustBeVar','''By'''));
                elseif ismember(by,var)
                    error(message('stats:fitrm:ByNotDistinct'));
                end
                var = [by var];
                iswithin = [bywithin iswithin];
            end

            % Get A and B matrices. One will represent a full factorial for
            % the selected factor; the other will be a default matrix.
            [A,Anames,Alevels] = emmMakeA(this,var,iswithin);
            [C,Cnames,Clevels] = emmMakeC(this,var,iswithin);
            if iswithin(end)
                numGroups = Clevels(end);
                [C,diffNames] = diffmatrix(C',Cnames,numGroups);
                C = C';
            else
                numGroups = Alevels(end);
                [A,diffNames] = diffmatrix(A,Anames,numGroups);
            end
            reorder = [];
            if isby
                if iswithin(end) && ~isempty(Anames)
                    % Differences involve Cnames; include Anames also
                    numA = size(Anames,1);
                    numDiff = size(diffNames,1);
                    rows = repmat((1:numA)',1,numDiff)';
                    diffNames = [Anames(rows(:),:),repmat(diffNames,numA,1)];
                elseif ~iswithin(end) && ~isempty(Cnames)
                    % Differences involve Anames; include Cnames also
                    numC = size(Cnames,1);
                    numDiff = size(diffNames,1);
                    rows = repmat((1:numC)',1,numDiff)';
                    diffNames = [Cnames(rows(:),:),repmat(diffNames,numC,1)];
                    numResults = size(diffNames,1);
                    if numC>0
                        reorder = reshape(1:numResults,numC,numResults/numC)';
                        reorder = reorder(:);
                    end
                end
            end

            % Carry out the emm calculation using this A and C
            [mn,se] = meanCalculate(this,A,C);
            mn = mn(:);
            se = se(:);
            t = mn./se;
            df = this.DFE;

            % Get critical values and adjusted p-values
            [crit,pval] = internal.stats.getcrit(ctype, alpha, df, numGroups, t);

            % Produce table of differences, standard errors, etc.
            width = se * crit;
            results = [mn,se,pval,mn-width,mn+width];
            if ~isempty(reorder)
                results = results(reorder,:);
            end
            tbl = [diffNames, ...
                   array2table(results,...
                         'VariableNames',{'Difference' 'StdErr' 'pValue' 'Lower' 'Upper'})];

        end
    end
    
    % --- Methods to support property access ---
    methods
        function this = set.WithinDesign(this,w)
        % This is a dependent property because its set method needs to
        % examine other properties including Y and ResponseNames.
            if ~isempty(w)
                ny = size(this.Y,2);
                oldnames = this.BetweenDesign.Properties.VariableNames;
                check = true;
                if isa(w,'dataset')
                    w = dataset2table(w);
                elseif ~isa(w,'table')
                    if isrow(w) && ny>1
                        w = w(:);
                    end
                    if iscolumn(w)
                        varnames = {'Time'};
                    else
                        varnames = internal.stats.numberedNames('w',1:size(w,2));
                    end
                    varnames = matlab.lang.makeUniqueStrings(varnames,oldnames);
                    if ~ismatrix(w)
                        error(message('stats:fitrm:WithinMatrix'));
                    end
                    w = array2table(w,'VariableNames',varnames);
                    check = false;
                end
                if size(w,1)~=ny
                    error(message('stats:fitrm:WithinBadLength',size(w,1),ny));
                end
                if check && any(ismember(w.Properties.VariableNames,oldnames))
                    error(message('stats:fitrm:WithinBadNames'));
                end
                w.Properties.RowNames = this.ResponseNames;
                if any(any(ismissing(w)))
                    error(message('stats:fitrm:NoMissing','WithinDesign'));
                end
            end
            this.WithinDesign_ = w;
        end
        function wd = get.WithinDesign(this)
            wd = this.WithinDesign_;
        end
        function this = set.WithinModel(this,w)
            if ~isempty(w) && ~(ischar(w) && ismember(w,{RepeatedMeasuresModel.SeparateMeans,RepeatedMeasuresModel.OrthogonalContrasts}))
                makeTestC(this,w); % just check for error
            end
            this.WithinModel = w;
        end
        function w = get.WithinFactorNames(this)
            if isempty(this.WithinDesign)
                w = {};
            else
                w = this.WithinDesign.Properties.VariableNames;
            end
        end
        function w = get.BetweenModel(this)
            w = this.Formula.LinearPredictor;
        end
        function w = get.BetweenFactorNames(this)
            factorCols = any(this.Terms,1);
            w = this.BetweenDesign.Properties.VariableNames(factorCols);
        end
        function w = get.ResponseNames(this)
            w = this.BetweenDesign.Properties.VariableNames(this.ResponseColumns);
        end
        function c = get.Covariance(this)
            rnames = this.ResponseNames;
            c = array2table(this.Cov,'VariableNames',rnames, 'RowNames',rnames);
            c.Properties.Description = getString(message('stats:fitrm:EstimatedCovariance'));
        end
        function c = get.Coefficients(this)
            rnames = this.ResponseNames;
            cnames = this.CoefNames;
            c = array2table(this.B,'VariableNames',rnames, 'RowNames',cnames);
            c.Properties.Description = getString(message('stats:fitrm:CoefficientEstimates'));
        end
        function c = get.DesignMatrix(this)
            c = this.X;
        end
    end
    
    % --- Static hidden method to support the FITRM function ---
    methods(Static, Hidden)
        function this = fit(ds,model,varargin) % fit method to be called by fitrm
%FIT Fit repeated measures model.
%   Do not call this function directly. Call FITRM instead.
%
%   See also FITRM.
            okargs =   {'WithinDesign' 'WithinModel'};
            defaults = {''             RepeatedMeasuresModel.SeparateMeans};
            [withindesign,withinmodel] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});
            
            if isa(ds,'dataset')
                ds = dataset2table(ds);
            elseif ~isa(ds,'table')
                error(message('stats:fitrm:TableRequired'));
            end

            % Get response columns from the formula. Make sure they are
            % valid.
            this = RepeatedMeasuresModel();
            varNames = ds.Properties.VariableNames;
            formula = classreg.regr.MultivariateLinearFormula(model,varNames);
            responseCols = names2cols(formula.ResponseName,varNames);
            f = @(var) isvector(var) && isnumeric(var) && isreal(var);
            ok = varfun(f,ds,'InputVariables',responseCols,'OutputFormat','uniform');
            if ~all(ok)
                responseCols = find(responseCols);
                bad = varNames(responseCols(~ok));
                error(message('stats:fitrm:BadResponse',bad{1}));
            end

            % Rows with missing values cannot be used in the fit
            missing = findMissing(ds,formula.Terms,responseCols);
            
            this.ResponseColumns = responseCols;
            Ymat = ds{~missing,responseCols};

            % Create a cell array to specify the valid values for
            % categorical predictors.
            iscat = varfun(@internal.stats.isDiscreteVar,ds,'OutputFormat','uniform');
            nvars = size(ds,2);
            vrange = cell(nvars,1);
            for i = 1:nvars
                vrange{i} = getVarRange(ds.(varNames{i}),iscat(i),missing);
            end

            % Compute the design matrix and related information from the
            % design matrix and the terms specified by the formula
            [Xmat,terms,~,coefTerms,coefNames,termNames] ...
                = classreg.regr.modelutils.designmatrix(ds(~missing,:),...
                'Model',formula.Terms, ...
                'DummyVarCoding','effects', ...
                'CategoricalVars',iscat, ...
                'CategoricalLevels',vrange);
             if rank(Xmat)<size(Xmat,2)
                 error(message('stats:fitrm:NotFullRank'));                        
             end
           
            if any(any(terms(:,responseCols)))
                factorCols = any(terms,1);
                bad = varNames(responseCols & factorCols);
                error(message('stats:fitrm:ResponseAndPredictor',bad{1}));
            end

            % Use the design matrix to carry out the fit, and compute
            % information needed to store in the object
            opt.RECT = true;
            Bmat = linsolve(Xmat,Ymat,opt);
            
            Resid = Ymat-Xmat*Bmat;
            dfe = size(Xmat,1)-size(Xmat,2);
            Covar = (Resid'*Resid)/dfe;
            
            % Term averages are required to computed expected marginal
            % means and to perform multiple comparisons of them
            this.TermAverages = calcTermAverages(ds,Xmat,formula.Terms,coefTerms,iscat);
            
            this.X = Xmat;
            this.Y = Ymat;
            this.Missing = missing;
            this.B = Bmat;
            this.Cov = Covar;
            this.CoefNames = coefNames;
            this.CoefTerms = coefTerms;
            this.TermNames = termNames;
            this.Formula = formula;
            this.BetweenDesign = ds;
            this.IsCat = iscat;
            this.VariableRange = vrange;
            this.DFE = dfe;
            this.Terms = terms;

            % Pre-compute Mauchly's test based on a full-rank test
            [this.Mauchly,this.Epsilon] = mauchlyTest(this.Cov,size(Xmat,1),rank(Xmat));
            
            if isempty(withindesign)
                withindesign = 1:size(Ymat,2);
            end
            this.WithinDesign = withindesign;
            this.WithinModel = withinmodel;
        end
    end
    
    % --- Protected methods ---
    methods(Access=protected)
        function this = RepeatedMeasuresModel() % constructor not to be called directly
        end
        function [mn,se] = meanCalculate(this,A,C) % calculate mean of ABC
            % Calculate estimated mean and standard error. This is used
            % both for expected marginal means and for predictions. Each
            % involves creating A and C matrices, then computing A*B*C with
            % its standard error.
            mn = (A*this.B*C)';

            % Compute the standard error if requested
            if nargout>=2
                S = this.Cov;
                XtX = this.X'*this.X;
                se = sqrt(sum(A.*(XtX\A')',2) * sum(C'.*(S*C)',2)')';
            end
        end
        function tbl = ranovastats(this,C,Acell,Anames,D,timename) % ranova statistics
            Beta = this.Coefficients{:,:};
            XX = this.DesignMatrix;
            [N,k] = size(XX);
            
            % Compute error matrix E
            % Would do [C,R]=qr(C,0) & convert D to D/R if D were not zero
            [C,~] = qr(C,0);
            E = (N-k) * C'*this.Covariance{:,:}*C;
            
            % Prepare variables for anova table
            r = size(C,2);
            nterms = numel(Acell);
            dfe = (N-k)*r;
            DF    = [ones(nterms,1); dfe];
            SumSq = [zeros(nterms,1); trace(E)];
            MeanSq = SumSq./DF;
            F = MeanSq/MeanSq(end);
            pValue = .5*ones(nterms+1,1);
            pValueGG = .5*ones(nterms+1,1);
            pValueHF = .5*ones(nterms+1,1);
            pValueLB = .5*ones(nterms+1,1);
            Eps = epsilon(this,C);
            for j=1:nterms
                A = Acell{j};
                s = size(A,1);
                
                % Compute hypothesis matrix H
                % H = (A*Beta*C - D)'*inv(A*inv(X'*X)*A')*(A*Beta*C - D);
                H = makeH(A,Beta,C,D,XX);
                
                % Compute test statistic.
                SumSq(j) = trace(H);
                DF(j) = s*r;
                MeanSq(j) = SumSq(j)/DF(j);
                F(j) = MeanSq(j)/MeanSq(end);
                [pValue(j),pValueGG(j),pValueHF(j),pValueLB(j)] ...
                    = pValueCorrections(F(j),DF(j),dfe,Eps);
            end
            absent = ((1:nterms+1)'==nterms+1);
            
            tbl = table(SumSq,DF,MeanSq);
            tbl.F = internal.stats.DoubleTableColumn(F,absent);
            tbl.pValue = internal.stats.DoubleTableColumn(pValue,absent);
            tbl.pValueGG = internal.stats.DoubleTableColumn(pValueGG,absent);
            tbl.pValueHF = internal.stats.DoubleTableColumn(pValueHF,absent);
            tbl.pValueLB = internal.stats.DoubleTableColumn(pValueLB,absent);
            
            if ~isequal(timename,'(Intercept)')
                Anames = strcat(Anames,':',timename);
                errorname = sprintf('%s(%s)','Error',timename);
            else
                errorname = 'Error';
            end
            tbl.Properties.RowNames = [Anames, errorname];
        end
        function [grp,iswithin] = getgroup(this,grp,catonly) % process grouping variable argument
            [tf,grp] = internal.stats.isStrings(grp);
            if ~tf
                error(message('stats:fitrm:GroupingNotCell'));
            end
            [tf,iswithin,isbetween,bloc] = isVariable(this,grp);
            if ~all(tf)
                bad = grp(~(iswithin | isbetween));
                error(message('stats:fitrm:GroupingNotRecognized',bad{1}));
            end
            if nargin>=3 && catonly
                iscat = this.IsCat;
                if any(~iscat(bloc(isbetween)))
                    error(message('stats:fitrm:GroupingNotCategorical'));
                end
            end
        end
        function [tf,iswithin,isbetween,bloc] = isVariable(this,var) % check for valid predictor variable
            iswithin = ismember(var,this.WithinFactorNames);
            [isbetween,bloc] = ismember(var,this.BetweenDesign.Properties.VariableNames);
            isresponse = ismember(var,this.ResponseNames);
            isbetween(isresponse) = false;
            tf = (iswithin | isbetween);
        end
        function [withindesign,time] = verifyWithinDesign(this,withinmodel,withindesign) % new design compatible with old
            time = [];
            oldw = this.WithinDesign;
            switch(withinmodel)
                case RepeatedMeasuresModel.OrthogonalContrasts
                    % Original design must be a single numeric vector
                    if size(oldw,2)~=1
                        error(message('stats:fitrm:NoTimeProperty'));
                    end
                    time = oldw{:,1};
                    if ~(isnumeric(time) && isvector(time) && isreal(time))
                        error(message('stats:fitrm:NoTimeProperty'));
                    end
                    
                    % New design must be compatible
                    if size(withindesign,2)~=1
                        error(message('stats:fitrm:NoTimeArgument'));
                    end
                    if istable(withindesign)
                        withindesign = withindesign{:,1};
                    end
                    if ~(isnumeric(withindesign) && isvector(withindesign) && isreal(withindesign))
                        error(message('stats:fitrm:NoTimeArgument'));
                    end
                    withindesign = withindesign(:);
                otherwise
                    [ok,withindesign] = verifyDesign(oldw,withindesign);
                    if ~ok
                        error(message('stats:fitrm:IncompatibleWithinDesign'));
                    end
            end
        end

        % --- Method to support customized display ---
        function group = getPropertyGroups(~) % properties and groups for display
            % Define three groups of properties, with titles
            titles = cell(3,1);
            titles{1} = sprintf('%s:',getString(message('stats:fitrm:DisplayBetween')));
            plist1 = {'BetweenDesign', 'ResponseNames', 'BetweenFactorNames', 'BetweenModel'};
            
            titles{2} = sprintf('%s:',getString(message('stats:fitrm:DisplayWithin')));
            plist2 = {'WithinDesign', 'WithinFactorNames', 'WithinModel'};
            
            titles{3} = sprintf('%s:',getString(message('stats:fitrm:DisplayEstimates')));
            plist3 = {'Coefficients', 'Covariance'};
            
            % Make group names bold if that is supported
            if matlab.internal.display.isHot
                for j=1:3
                    titles{j} = sprintf('<strong>%s</strong>',titles{j});
                end
            end
            
            group(1) = matlab.mixin.util.PropertyGroup(plist1,titles{1});
            group(2) = matlab.mixin.util.PropertyGroup(plist2,titles{2});
            group(3) = matlab.mixin.util.PropertyGroup(plist3,titles{3});
        end
    end
end

% ---------- local functions -----------
function missing = findMissing(t,terms,responseCols) % find rows with missing values

useCols = any(terms,1);
useCols(responseCols) = true;
missing = any(ismissing(t(:,useCols)),2);
end

function ds = emmMakeDataset(Bff,Wff,grp,m,se,dfe,alpha) % make dataset of emm values
% Create table of means and standard errors labeled by the Between and
% Within factors
if isempty(Bff)
    ds = Wff;
elseif isempty(Wff)
    ds = Bff;
else
    nb = size(Bff,1);
    nw = size(Wff,1);
    Brows = kron((1:nb)',ones(nw,1));
    Wrows = repmat((1:nw)',nb,1);
    ds = [Bff(Brows,:),Wff(Wrows,:)];
end
ds = ds(:,grp);

t = -tinv(alpha/2,dfe);
ds.Mean = m;
ds.StdErr = se;
ds.Lower = m - t*se;
ds.Upper = m + t*se;

ds = sortrows(ds,grp);
end

function [X,termnames] = getTermAverage(this,termNums) % compute term averages
iscat = this.IsCat;
savedAverages = this.TermAverages;

% Find continuous part of each term among the saved averages
terms = this.Formula.Terms;
terms(:,iscat) = 0;
[~,loc] = ismember(terms,savedAverages{1},'rows');

% Expand to full term size, then to full design matrix size
termAverages = savedAverages{2}(loc)';

coefTerms = this.CoefTerms;
coefTerms = coefTerms(ismember(coefTerms,termNums));
X = termAverages(coefTerms);

if nargout>=2
    termnames = classreg.regr.modelutils.terms2names(savedAverages{1},this.Formula.VariableNames);
    termnames = termnames(loc(coefTerms));
end
end

function cols = names2cols(str,varNames) % get response columns from names
ranges = textscan(str,'%s','delimiter',' ,','MultipleDelimsAsOne',true);
ranges = ranges{1};
cols = false(1,length(varNames));
for j=1:length(ranges)
    range = ranges{j};
    loc = find(range=='-');
    if isempty(loc)
        cols(strcmp(range,varNames))= true;
    else
        loc1 = find(strcmp(range(1:loc-1),varNames));
        loc2 = find(strcmp(range(loc+1:end),varNames));
        cols(min(loc1,loc2):max(loc1,loc2)) = true;
    end
end
end

function range = getVarRange(v,asCat,excl) % variable range for design matrix calculations
v(excl,:) = [];
if asCat % create a list of the unique values
    if isa(v,'categorical')
        % For categorical classes, get the values actually present in the
        % data, not the set of possible values.
        range = unique(v(:));
        range = range(~isundefined(range)); % empty if NaNs
        
    else
        % For classes other than categorical, the list of unique values is
        % also the list of possible values.  But get it in the same order as
        % what grp2idx defines for each class.
        [~,~,range] = grp2idx(v); % leaves NaN or '' out of glevels
        
    end
    if ~ischar(range)
        range = range(:)'; % force a row
    end
elseif isnumeric(v) || islogical(v) % find min and max
    range = [min(v,[],1)  max(v,[],1)]; % ignores NaNs unless all NaNs
else
    range = NaN(1,2);
end
end

function ds = manovastats(X,A,B,C,D,SSE,varargin) % all manova tests
% Compute four statistics for testing
%     A*B*C = D
% where B is the coefficient matrix and the other matrices are inputs.
%
% The varargin parameter consists of the between-subject names and the
% within-subject names.
if nargin<5
    D = 0; % default is to test for zero
end

if ~iscell(A)
    A = {A};
end
if ~iscell(C)
    C = {C};
end
na = numel(A);
nc = numel(C);

dscell = cell(na,nc);
for k = 1:na
    for j=1:nc
        dscell{k,j} = fourstats(X,A{k},B,C{j},D,SSE);
    end
end
if na*nc>1
    ds = combinemanovatables(dscell,varargin{:});
else
    ds = dscell{1};
end
end

function ds = fourstats(X,A,B,C,D,SSE) % four manova statistics for one test

% Hypothesis matrix H
% H = (A*Beta*C - D)'*inv(A*inv(X'*X)*A')*(A*Beta*C - D);
% q = rank(Z);
[H,q] = makeH(A,B,C,D,X);

% Error matrix E
E = C'*SSE*C;

checkHE(H,E);

p = rank(E+H);
s = min(p,q);
v = size(X,1)-rank(X);
if p^2+q^2>5
    t = sqrt( (p^2*q^2-4) / (p^2+q^2-5));
else
    t = 1;
end
u = (p*q-2)/4;
r = v - (p-q+1)/2;
m = (abs(p-q)-1)/2;
n = (v-p-1)/2;

% ~~~ Wilks' Lambda = L
% Formally, L = |E| / |H+E|, but it is more convenient to compute it using
% the eigenvalues from a generalized eigenvalue problem
lam = eig(H,E);
mask = (lam<0) & (lam>-100*eps(max(abs(lam))));
lam(mask) = 0;
L_df1 = p*q;
L_df2 = r*t-2*u;
if isreal(lam) && all(lam>=0) && L_df2>0
    L = prod(1./(1+lam));
else
    L = NaN;
    L_df2 = max(0,L_df2);
end
L1 = L^(1/t);
L_F = ((1-L1) / L1) * (r*t-2*u)/(p*q);
L_rsq = 1-L1;

% ~~~ Pillai's trace = V
% Formally, V = trace( H * (H+E)^-1 ), but again it can be reduced to a
% generalized eigenvalue problem
theta = eig(H,H+E);
V = sum(theta);
if s>V
    V_F = ( (2*n+s+1)/(2*m+s+1)) * V / (s - V);
    V_rsq = V/s;
else
    V_F = NaN;
    V_rsq = NaN;
end
V_df1 = s*(2*m+s+1);
V_df2 = s*(2*n+s+1);

% ~~~ Hotelling-Lawley trace = U
% U = trace( H * E^-1 ) but it can also be written as the sum of the
% eigenvalues that we already obtained above
if isreal(lam) && all(lam>=0)
    U = sum(lam);
else
    U = NaN;
    n(n<0) = NaN;
end
b = (p+2*n)*(q+2*n) / (2*(2*n+1)*(n-1));
c = (2 + (p*q+2)/(b-1))/(2*n);
if n>0
    U_F = (U/c) * (4+(p*q+2)/(b-1)) / (p*q);
else
    U_F = U * 2 * (s*n+1) / (s^2 * (2*m+s+1));
end
U_rsq = U / (U+s);
U_df1 = s*(2*m+s+1);
U_df2 = 2*(s*n+1);

% ~~~ Roy's maximum root = R
% This is the maximum eigenvalue already found, though sometimes the
% quantity defined as R_rsq below is used as the definition of Roy's
% statistic
if isempty(lam)
    R = NaN;
else
    R = max(lam);
end
r = max(p,q);
R_F = R * (v-r+q)/r;
R_rsq = R/(1+R);
R_df1 = r;
R_df2 = v-r+q;

% Create a table of all four statistics with their corresponding R-square
% and F values
Statistic = categorical({'Pillai' 'Wilks' 'Hotelling' 'Roy'}');
Value = [V;L;U;R];
F = [V_F;L_F;U_F;R_F];
RSquare = [V_rsq; L_rsq; U_rsq; R_rsq];
df1 = [V_df1;L_df1;U_df1;R_df1];
df2 = [V_df2;L_df2;U_df2;R_df2];
ds = table(Statistic, Value, F, RSquare,df1,df2);
ds.pValue = fcdf(F, df1, df2, 'upper');
end

function [tbl,B] = contrastanova(d,formula,Y,Ynames,drop,yname) % anova calculations on multiple inputs
if nargin<5
    drop = false;
end
if nargin<6
    yname = 'y';
end

ny = size(Y,2);
C = cell(1,ny);
B = zeros(0,ny);
for j=1:ny
    d.(yname) = Y(:,j);
    lm = LinearModel.fit(table2dataset(d),formula,'DummyVarCoding','effects');
    C{j} = anovawithconstant(lm);
    b = lm.Coefficients.Estimate;
    B(1:length(b),j) = b;
end
if drop
    C(1) = [];
    Ynames(1) = [];
end
tbl = combinemanovatables(C,{},Ynames);
end

function a = anovawithconstant(lm) % anova calculations on one input

% Get regular anova
a = anova(lm,'components',3);
if isa(a,'dataset')
    a = dataset2table(a);
end

% Reconstruct information for constant term
p0 = lm.Coefficients.pValue(1);
f0 = lm.Coefficients.tStat(1)^2;
mse = a.MeanSq(end);
ms0 = f0*mse;
ss0 = ms0;

% Append to anova table
a = a([1:end,end],:);
a.Properties.RowNames{end} = 'constant';
a.SumSq(end) = ss0;
a.DF(end) = 1;
a.MeanSq(end) = ms0;
a.F(end) = f0;
a.pValue(end) = p0;

% Put this one first
a = a([end,1:end-1],:);

% Reset the entries not to display
n = size(a,1);
absent = (1:n)' == n;
a.F = internal.stats.DoubleTableColumn(a.F,absent);
a.pValue = internal.stats.DoubleTableColumn(a.pValue,absent);
end

function D = combinemanovatables(C,bnames,wnames) % combine manova/anova tables
%COMBINEMANOVATABLES Combine manova or anova tables into a single table
%   D = combinemanovatables(C,b,w)
%         C = B-by-1 cell array of separate tables
%         b = cell array of B strings
%         w = cell array of W strings

% Measure inputs. Set bnames to be the row names of one table if no input
% value is given.
first = C{1};
nr = size(first,1);
if nargin<2 || isempty(bnames)
    bnames = first.Properties.RowNames;
    nr = 1;
    nb = size(first,1);
else
    nb = length(bnames);
end
if nargin<3 || isempty(wnames)
    wnames = {'Within'};
end
nw = length(wnames);

% Combine tables
B = numel(C);
D2 = repmat(C{1},nw,1);
numrows = size(C{1},1);
for k=2:B
    base = numrows*(k-1);
    D2(base+(1:numrows),:) = C{k};
end

% Create names for the Within column
wnames = wnames(:);
t = repmat(1:numel(wnames),nr*nb,1);
t = t(:);
Within = categorical(wnames(t));
D1 = table(Within);

% Create names for the Between column
Between = repmat(bnames(:)',nr,nw);
Between = categorical(Between(:));
D1.Between = Between;

% Combine names with data
if ismember('Within',D2.Properties.VariableNames)
    D2.Within = [];
end
D = [D1 D2];

% Remove observation names from result
D.Properties.RowNames = {};

% See if we need to suppress any entries
classes = varfun(@(v){class(v)},first,'OutputFormat','uniform');
dtccol = find(strcmp('internal.stats.DoubleTableColumn',classes));
if ~isempty(dtccol)
    vnames = first.Properties.VariableNames;
    absent = first.(vnames{dtccol(1)}).absent;
    absent = repmat(absent,nr*nw,1);
    for j=1:length(dtccol)
        vname = D.Properties.VariableNames{dtccol(j)+2};
        D.(vname) = internal.stats.DoubleTableColumn(D.(vname),absent);
    end
end
end

function [str,epsilon] = mauchlyTest(S,n,rx,C) % Mauchly's test of sphericity
% Mauchly's test for sphericity

p = size(S,1);
d = p-1;

% Orthonormal contrast matrix C
if nargin<4
    C = triu(ones(p,d));
    C(2:p+1:end) = -(1:d);
    s = sqrt(2*cumsum(1:d));
    C = bsxfun(@rdivide,C,s);
else
    [C,~] = qr(C,0);
end
d = size(C,2);

% Mauchly's W statistic
% W = det(T)/(trace(T)/p)^d     where T is C'*S*C and 
%                                     d=p-1 is the number of columns of T
%   = prod(lam) / mean(lam)^d   where lam is the eigenvalues of T
%   = prod( lam/mean(lam) )
lam = eig(C'*S*C);
lam = max(0,real(lam));
avglam = sum(lam)/d;
W = prod(lam/avglam);

% Epsilon correction factors
if nargout>=2
    Uncorrected = 1;
    if isempty(lam)
        GreenhouseGeisser = 1;
        HuynhFeldt = 1;
        LowerBound = 1;
    else
        LowerBound = 1/d;
        GreenhouseGeisser = min(1,max(LowerBound,...
            sum(lam)^2 / (d*sum(lam.^2))));
        
        % Huynh-Feldt correction
        % Original:
        %    HuynhFeldt = min(1, (n*d*GreenhouseGeisser-2) / (d*(n-rx)-d^2*GreenhouseGeisser))
        % Lecoutre modification:
        HuynhFeldt = min(1,max(LowerBound,...
            ((n-rx+1)*d*GreenhouseGeisser-2) / (d*(n-rx) - d^2*GreenhouseGeisser)));
    end
    epsilon = table(Uncorrected,GreenhouseGeisser,HuynhFeldt,LowerBound);
end

% Chi-square test

nr = n - rx;
dd = 1 - (2*d^2+d+2) / (6*d*nr);
ChiStat = -log(W) * dd * nr;
DF = max(0,d*(d+1)/2 - 1);
p1 = chi2cdf(ChiStat,DF,'upper');

pValue = p1;                        % first-order approximation

% Below is the second-order approximation recommended in T.W. Anderson,
% "An Introduction to Multivariate Statistical Analysis," p. 436 of the 3rd
% edition.
% p2 = chi2cdf(ChiStat,DF+4,'upper');
% w2 = (d+2).*(d-1).*(d-2).*(2*d.^3+6*d.^2+3*d+2)  ./  (288*d.^2.*nr^2.*dd.^2);
% pValue = p1 + w2.*(p2-p1);        % second-order approximation

str = table(W, ChiStat, DF, pValue);
end

function [pUnc,pGG,pHF,pLB] = pValueCorrections(F,df,dfe,Epsilon) % compute corrected p-values

% Compute uncorrected p-value
pUnc = fcdf(F, df, dfe, 'upper');

% Compute corrected values
e = Epsilon.GreenhouseGeisser(1);
pGG = fcdf(F, e*df, e*dfe, 'upper');

e = Epsilon.HuynhFeldt(1);
pHF = fcdf(F, e*df, e*dfe, 'upper');

e = Epsilon.LowerBound(1);
pLB = fcdf(F, e*df, e*dfe, 'upper');
end

function [ngroups,grpidx,grpname] = makegroup(group,this) % get grouping indices defined by grouping variable
if isempty(group)
    ngroups = 1;
    nsubjects = size(this.X,1);
    grpidx = ones(nsubjects,1);
    grpname = {};
    return
end

[tf,gvars] = internal.stats.isStrings(group);
if ~tf
    error(message('stats:fitrm:GroupingNotCell'))
end
if ~all(ismember(gvars,this.BetweenDesign.Properties.VariableNames))
    error(message('stats:fitrm:GroupingNotBetween'))
end
groups = cell(1,length(gvars));
for j = 1:length(gvars)
    gj = this.BetweenDesign.(gvars{j});
    groups{j} = gj(~this.Missing,:);
end
[grpidx,~,grpvals] = internal.stats.mgrp2idx(groups);
ngroups = size(grpvals,1);
grpname = cell(ngroups,1);
for k = 1:ngroups
    grpname{k} = sprintf('%s=%s',gvars{1},grpvals{k,1});
end
for j = 2:length(gvars)
    for k = 1:ngroups
        grpname{k} = sprintf('%s, %s=%s',grpname{k},gvars{j},grpvals{k,j});
    end
end
end

function c = calcTermAverages(ds,Xmat,terms,termCols,iscat) % calculate term averages

% Find terms representing continuous predictors only
contterms = terms;
contterms(:,iscat) = 0;
contterms = unique(contterms,'rows');

[tf,loc] = ismember(contterms,terms,'rows');
averages = zeros(size(contterms,1),1);

% Ordinarily we will already have these terms in the model, but in case we
% do not, compute the values of the missing terms
if any(~tf)
    Zmat = classreg.regr.modelutils.designmatrix(ds,'Model',contterms(~tf,:));
end

% Compute term averages from either X or Z
Zcol= 0;
for j=1:numel(tf)
    if tf(j)
        termloc = loc(j);
        avg = mean(Xmat(:,termloc==termCols));
    else
        Zcol = Zcol + 1;
        avg = mean(Zmat(:,Zcol));
    end
    averages(j) = avg;
end
c = {contterms averages};
end

function [ff,nlevels] = makeFullFactorial(design) % create full factorial
ngroups = size(design,2);
nlevels = zeros(1,ngroups);
levels = cell(1,ngroups);
grp = design.Properties.VariableNames;
for j=1:ngroups
    gj = design.(grp{j});
    if ischar(gj)
        gj = cellstr(gj);
    end
    u = unique(gj);
    levels{j} = u;
    nlevels(j) = size(u,1);
end
if ngroups>0
    indexdesign = fliplr(fullfact(fliplr(nlevels))); % last col varies fastest
end
ff = table();
for j=1:ngroups
    ff.(grp{j}) = levels{j}(indexdesign(:,j),:);
end
end

function [C,Cnames] = makeTestC(this,model,celloutput,dsin,numok) % C matrix for ABC=D test
%makeTestC Compute contrast matrix for testing model terms.
%    C = makeTestC(RM,MODEL,CELLOUTPUT,DSIN,NUMOK) returns a contrast matrix
%    C for testing the terms in the model MODEL for the within-subjects
%    factors. If CELLOUTPUT is TRUE, the output C is a cell array with a
%    separate entry for each term. If MODEL is 'separatemeans', the output
%    C is a matrix with contrasts for a test of equality for the separate
%    means. DSIN is the table to use in computing the contrasts, and it is
%    RM.WithinDesign by default. NUMOK is true if it is acceptable to
%    supply C as the model.
%
%    [C,CNAMES] = makeTestC(...) also returns a cell array CNAMES of names
%    for the contrast matrix or matrices.
if nargin<3
    celloutput = true;
end
if nargin<4 || isempty(dsin)
    dsin = this.WithinDesign;
end
if nargin<5
    numok = false;
end
ny = size(this.Y,2);
if numok && isnumeric(model)
    C = model;
    checkMatrix('WithinModel',C,ny,[]);
    Cnames = {getString(message('stats:fitrm:SpecifiedContrast'))};
    celloutput = false;    % just one "term"
elseif isequal(model,RepeatedMeasuresModel.SeparateMeans)
    % just test for equal means
    C = eye(ny-1,ny);
    C(ny:ny:end) = -1;
    C = C';
    Cnames = {'Constant'};
    celloutput = false;    % just one "term"
elseif isequal(model,RepeatedMeasuresModel.MeanResponse)
    C = ones(ny,1)/sqrt(ny);
    Cnames = {'Constant'};
    celloutput = false;    % just one "term"
elseif isequal(model,RepeatedMeasuresModel.OrthogonalContrasts)
    W = dsin;
    if size(W,2)~=1 || ~varfun(@isnumeric,W,'OutputFormat','uniform')
        error(message('stats:fitrm:NotNumericFactor'));
    end
    timename = W.Properties.VariableNames{1};
    time = W.(timename);
    W = fliplr(vander(time(:)));
    [C,~] = qr(W);
    Cnames = cell(1,ny);
    Cnames{1} = 'Constant';
    Cnames{2} = timename;
    for j=2:ny-1
        Cnames{j+1} = sprintf('%s^%d',timename,j);
    end
    
else
    ds = this.WithinDesign;
    yname = genvarname('y',ds.Properties.VariableNames);
    vnames = [ds.Properties.VariableNames {yname}];
    if ischar(model)
        try
            formula = classreg.regr.LinearFormula([yname ' ~ ' model],vnames);
        catch ME
            ME = addCause(MException(message('stats:fitrm:BadWithinModel')),ME);
            throw(ME);
        end
    else
        error(message('stats:fitrm:BadWithinModel'));
    end
        
    varNames = ds.Properties.VariableNames;
    if ~all(ismember(varNames,dsin.Properties.VariableNames))
        error(message('stats:fitrm:BadWithinModel'));
    else
        dsin = dsin(:,varNames);
    end
    iscat = varfun(@internal.stats.isDiscreteVar,ds,'OutputFormat','uniform');
    nvars = size(ds,2);
    vrange = cell(nvars,1);
    excl = [];
    for i = 1:nvars
        vrange{i} = getVarRange(ds.(varNames{i}),iscat(i),excl);
    end
    
    terms = formula.Terms(:,1:end-1);
    [C,~,~,termcols,Cnames,TermNames] = classreg.regr.modelutils.designmatrix(dsin,...
        'Model',terms, ...
        'DummyVarCoding','effects', ...
        'CategoricalVars',iscat, ...
        'CategoricalLevels',vrange);
end
if celloutput
    termsize = accumarray(termcols',1);
    C = mat2cell(C,size(C,1),termsize);
    Cnames = TermNames;
end
end

function [A,bnames] = makeTestA(this) % A matrix for ABC=D test
ncoefs = length(this.CoefTerms);
termsize = accumarray(this.CoefTerms',1);
A = mat2cell(eye(ncoefs),termsize);
bnames = this.TermNames;
end

function [A,bnames] = makeTestABy(this,by) % A matrix for ABC=D best by another vairable
if ~internal.stats.isString(by) || ~ismember(by,this.BetweenFactorNames)
    error(message('stats:fitrm:ByNotBetween'))
end

% use between-subjects model terms
byvar = this.BetweenDesign.(by);
byvar = byvar(~this.Missing,:);
if ischar(byvar)
    byvar = cellstr(byvar);
end
[~,vnum] = ismember(by,this.BetweenDesign.Properties.VariableNames);
if this.IsCat(vnum)
    u = this.VariableRange{vnum};
    if ischar(u)
        u = cellstr(u);
    end
else
    u = unique(byvar);
end
A = cell(length(u),1);
bnames = cell(length(u),1);
for j=1:length(u)
    if iscell(u)
        uj = u{j};
        rows = strcmp(uj,byvar);
    else
        uj = u(j);
        rows = byvar==uj;
    end
    A{j} = mean(this.X(rows,:),1);
    if isnumeric(uj) || islogical(uj)
        uj = num2str(uj);
    else
        uj = char(uj);
    end
    bnames{j} = sprintf('%s=%s',by,uj);
end
end

function [Xnew,Bff,nlevels,cap] = emmMakeA(this,grp,iswithin) % A matrix for emm calculation
% Create a full factorial design for the Between variables
[~,vnums] = ismember(grp(~iswithin),this.BetweenDesign.Properties.VariableNames);
[~,order] = sort(vnums);
[Bff,nlevels] = makeFullFactorial(this.BetweenDesign(~this.Missing,grp(~iswithin)));
Bvars = ismember(this.BetweenDesign.Properties.VariableNames,grp);

% Find other between-subject categorical variables
iscat = this.IsCat;
othercat = find(iscat);
bnames = this.BetweenDesign.Properties.VariableNames;
othercat(ismember(bnames(othercat),grp)) = [];

% Determine which terms involve the unused categorical variables. Since we
% have used effects coding, we can omit these terms.
terms = this.Formula.Terms;
temp = any(terms(:,othercat),2)';
keepTerms = find(~temp);

% For each remaining term, compute the design matrix columns for the
% categorical part of the term.
catpart = terms(keepTerms,:);
catpart(:,~iscat) = 0;
if isempty(Bff)
    X1mat = 1;
else
    X1mat = classreg.regr.modelutils.designmatrix(Bff(:,order),...
        'Model',catpart(:,Bvars), ...
        'DummyVarCoding','effects', ...
        'CategoricalVars',iscat(Bvars), ...
        'CategoricalLevels',this.VariableRange(Bvars));
end

% Multiply each by the average over the continuous part of the term.
[X2row,avgNames] = getTermAverage(this,keepTerms);
Xnew = zeros(size(X1mat,1),size(this.X,2));
Xnew(:,ismember(this.CoefTerms,keepTerms)) = bsxfun(@times,X1mat,X2row);

% Create caption with covariate values if necessary.
covariates = ~strcmp(avgNames,'(Intercept)');
if any(covariates)
    cap = strcat(avgNames(covariates),'=',num2str(X2row(covariates)','%-g'));
    cap = sprintf(', %s',cap{:});
    cap = getString(message('stats:fitrm:MeansComputedWith',cap(3:end)));
else
    cap = '';
end
end

function [M,Wff,nlevels] = emmMakeC(this,grp,iswithin) % C matrix for emm calculation
if ~any(iswithin)
    ny = size(this.B,2);
    M = ones(ny,1)/ny;
    Wff = [];
    nlevels = 1;
else
    [Wff,nlevels] = makeFullFactorial(this.WithinDesign(:,grp(iswithin)));
    w = this.WithinDesign;
    cols = ismember(this.WithinFactorNames,grp);
    [tf,loc] = ismember(w(:,cols),Wff);
    if any(~tf)
        error(message('stats:fitrm:CombinationsMissing'));
    end
    M = zeros(size(w,1),size(Wff,1));
    for j = 1:size(Wff,1)
        t = (loc==j);
        M(:,j) = t/sum(t);
    end
end
end

function [terms,iscat,vrange] = dataset2terms(this,ds) % compute terms from data

% Get term and variable info from the model
rmvars = this.BetweenDesign.Properties.VariableNames;
terms = this.Formula.Terms;
iscat = this.IsCat;
vrange = this.VariableRange;

% Get variable names from the dataset; done if they match
dsvars = ds.Properties.VariableNames;
if isequal(dsvars,rmvars)
    if ~verifyDesign(this.BetweenDesign,ds);
        error(message('stats:fitrm:IncompatibleBetweenDesign'));
    end
    return
end

% Make sure all required variables are present and of the right type
bfvars = this.BetweenFactorNames;
[ok,dsidx] = ismember(bfvars,dsvars);
if ~all(ok)
    vname = bfvars(~ok);
    error(message('stats:fitrm:BetweenFactorMissing',vname{1}));
end
[~,rmidx] = ismember(bfvars,rmvars);
nvars = size(ds,2);
if ~verifyDesign(this.BetweenDesign(:,rmidx),ds(:,dsidx));
    error(message('stats:fitrm:IncompatibleBetweenDesign'));
end

% Modify terms etc. to match input dataset
newterms = zeros(size(terms,1),nvars);
newterms(:,dsidx) = terms(:,rmidx);
terms = newterms;

newcat = false(nvars,1);
newcat(dsidx) = iscat(rmidx);
iscat = newcat;

newrange = repmat(vrange(1),nvars,1);
newrange(dsidx) = vrange(rmidx);
vrange = newrange;
end

function D = alldiff(n) % contrast matrix for all differences
% Compute contrast matrix for all differences among N elements, including
% differences in both directions

% Create a matrix of (I,J) pairs with I~=J
D = fliplr(fullfact([n n]));  % all (I,J) pairs
D(D(:,1)==D(:,2),:) = [];     % remove pairs with I=J
end

function [D,Dnames] = diffmatrix(A,Anames,numGroups) % create matrix to represent differences
% Create matrix D representing differences between all ordered pairs of
% distinct rows of A, or each set of numGroups rows of A if we're doing
% this 'by' one or more grouping variables
if nargin<3
    numGroups = size(A,1);
end
numBy = size(A,1)/numGroups;
pairs = alldiff(numGroups);
one = pairs(:,1);
two = pairs(:,2);
if numBy>1
    one = bsxfun(@plus,one,numGroups*(0:numBy-1));
    one = one(:);
    two = bsxfun(@plus,two,numGroups*(0:numBy-1));
    two = two(:);
end
D = A(one,:) - A(two,:);

varnames = Anames.Properties.VariableNames{end};
first = Anames(one,end);
first.Properties.VariableNames = {[varnames,'_1']};
second = Anames(two,end);
second.Properties.VariableNames = {[varnames,'_2']};
Dnames = [Anames(one,1:end-1),first, second];
end

function [cmap,markers,styles] = regularizePlotArgs(cmap,markers,styles,ngroups) % convert plot args to standard form
if nargin<4
    ngroups = 1;
end
if isempty(cmap)
    cmap = lines(ngroups);
else
    cmap = internal.stats.colorStringToRGB(cmap);
end
if ischar(markers)
    markers = {markers};
end
if ischar(styles)
    styles = {styles};
end
end

function checkAlpha(alpha) % check alpha parameter
if ~isscalar(alpha) || ~isnumeric(alpha) || ~isreal(alpha) ...
                    || ~isfinite(alpha)  || alpha<=0 || alpha>=1
    throwAsCaller(MException(message('stats:fitrm:BadAlpha')))
end
end

function checkMatrix(name,A,rows,cols,okscalar) % check for value coeftest matrix
% Check that the "name" matrix A is the requested number of rows and
% columns. The optional fourth argument indicates whether, if both rows and
% columns are specified that a scalar is also accepted.
msg = [];

if any(any(isnan(A)))
    msg = message('stats:fitrm:NoMissing',name);
end

checkA = isempty(rows);
checkC = isempty(cols);

% The following enforce the specific constraints needed in this file:
if checkA                     % checking the A matrix
    if ~isnumeric(A) || ~ismatrix(A) || size(A,2)~=cols
        msg = message('stats:fitrm:MatrixWithCols',name,cols);
    elseif any(all(A==0,2))
        msg = message('stats:fitrm:BadAMatrix');
    end
elseif checkC                 % checking the C matrix
    if ~isnumeric(A) || ~ismatrix(A) || size(A,1)~=rows
        msg = message('stats:fitrm:MatrixWithRows',name,rows);
    end
else                          % checking the D matrix or similar
    if ~isnumeric(A) || ~(isscalar(A) || isequal(size(A),[rows,cols]))
        if nargin<5 || okscalar
            msg = message('stats:fitrm:MatrixWithSize',name,rows,cols);
        else
            msg = message('stats:fitrm:MatrixWithSizeNoScalar',name,rows,cols);
        end
    end
end
if ~isempty(msg)
    throwAsCaller(MException(msg));
end

if (checkA||checkC) && rank(A) < min(size(A))
    error(message('stats:fitrm:RankDefA',name));
end
end

function [H,q] = makeH(A,B,C,D,X) % Make hypothesis matrix H
% H = (A*Beta*C - D)'*inv(A*inv(X'*X)*A')*(A*Beta*C - D);
d = A*B*C - D;
[~,RX] = qr(X,0);
XA = A/RX;
Z = XA*XA';
H = d'*(Z\d);  % note Z is often scalar or at least well-conditioned

if nargout>=2
    q = rank(Z);
end
end

function checkHE(H,E) % Check that Error and Hypothesis matrices are okay
msg = [];
if ~all(all(isfinite(H)))
    msg = message('stats:fitrm:BadHMatrix');
elseif ~all(all(isfinite(E)))
    msg = message('stats:fitrm:BadEMatrix');
end
if ~isempty(msg)
    throwAsCaller(MException(msg));
end
end

function [ok,new] = verifyDesign(old,new)
oldnum = varfun(@(x)isnumeric(x),old,'OutputFormat','uniform');
ok = false;
if ~istable(new)
    % A matrix is okay for a pure numeric design
    if size(new,2)~=length(oldnum) || ~all(oldnum) || ~ismatrix(new) || ~isnumeric(new)
        ok = false;
        return
    end
    new = array2table(new,'VariableNames',old.Properties.VariableNames);
else
    % Tables must have the same length and variable classes
    oldclass = varfun(@(v){class(v)},old,'OutputFormat','uniform');
    newclass = varfun(@(v){class(v)},new,'OutputFormat','uniform');
    if isequal(oldclass,newclass)
        ok = true;
        return
    end
    if numel(oldclass)~=numel(newclass)
        return
    end
    
    % Allow char input for categorical classes
    for j=1:numel(oldclass)
        if isequal(oldclass{j},newclass{j})
            continue
        end
        vold = old{:,j};
        if isa(vold,'categorical')
            vnew = new{:,j};
            if    isa(vnew,'categorical') ...
                    || (ismatrix(vnew) && ischar(vnew)) ...
                    || iscellstr(vnew)
                continue
            end
        end
        return;
    end
end
ok = true;
return
end
