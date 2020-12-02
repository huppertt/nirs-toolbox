function s = fitrm(ds,model,Noise,pmax,varargin)
%FITRM Fit repeated measures model.
%   RM = FITRM(T,MODELSPEC) fits the model specified by MODELSPEC to
%   data in the table T, and returns the repeated measures model RM.
%   T is a table containing the values of the response variables and
%   the between-subject factors to be used as predictors in the model.
%
%   MODELSPEC specifies the response variable names and model terms
%   as a string such as 'y1-y5 ~ x1 + x2 + x3*x4'. The string contains
%   a response specification, followed by '~', followed by one or more
%   terms joined by '+' or '-'. The following are rules for specifying
%   responses:
%            Y1,Y2,Y3  the specified list of variables
%            Y1-Y5     all table variables from Y1 to Y5
%
%   The following are rules for specifying terms:
%            A + B    term A and term B
%            A - B    term A but without term B
%            A:B      the product of A and B
%            A*B      A + B + A:B
%            A^2      A + A:A
%            ()       grouping of terms
%
%   Variables used in model terms are treated as categorical if they are
%   categorical (nominal or ordinal), logical, char arrays, or cell
%   arrays of strings.
%       
%   RM = FITRM(T,MODELSPEC,NAME1,VAL1,...) also specifies one or more
%   of the following name/value pairs:
% 
%        'WithinDesign' A design for within-subject factors, specified
%                       as one of the following:
%               V  - a numeric vector of length R, where R is the
%                    number of repeated measures, specifying the R
%                    measurement times. 
%               M  - an R-by-K numeric matrix of the values of K
%                    within-subject factors w1,w2,...,wK. All K variables
%                    are treated as continuous.
%               TW - a table with R rows and K variables containing
%                    the values of K within-subject factors.
%
%        'WithinModel'  A model WMODEL specifying the within-subjects
%                       hypothesis test. WMODEL may be any of the following:
%
%            'separatemeans'
%                    Compute a separate mean for each group
%            'orthogonalcontrasts'
%                    Valid only when then the within-factor design
%                    consists of a single numeric factor. If that factor
%                    is Time, specifies a model consisting of orthogonal
%                    polynomials up to order Time.^(R-1).
%            WMODELSPEC
%                    A text string that specifies the terms involving the
%                    within-subject factors, following the rules for the
%                    part of the MODELSPEC to the right of the '~'.
%
%   Example: Model the repeated measurements y1-y8 using a between-subject
%            design that includes one interaction term plus additional
%            main effects, plus a within-subject design with two factors.
%      load repeatedmeas
%      R = fitrm(between, 'y1-y8 ~ Group*Gender+Age+IQ', 'WithinDesign',within)
%
%   See also RepeatedMeasuresModel.

%   Copyright 2013-2017 The MathWorks, Inc.

internal.stats.checkNotTall(upper(mfilename),0,ds,model,varargin{:});

if nargin > 1
    model = convertStringsToChars(model);
end

if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

s = nirs.math.RepeatedMeasuresModel.fit(ds,model,Noise,pmax,varargin{:});


return