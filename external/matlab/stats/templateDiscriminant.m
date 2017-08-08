function temp = templateDiscriminant(varargin)
%TEMPLATEDISCRIMINANT Create a discriminant template.
%   T=TEMPLATEDISCRIMINANT() returns a discriminant template suitable for
%   use in the FITENSEMBLE and FITCECOC functions.
%
%   T=TEMPLATEDISCRIMINANT('PARAM1',val1,'PARAM2',val2,...) specifies
%   optional parameter name/value pairs:
%
%       'DiscrimType'           - A case-insensitive string with the type
%                                 of the discriminant analysis. Specify as
%                                 one of: 'linear', 'pseudolinear',
%                                 'diaglinear', 'quadratic',
%                                 'pseudoquadratic' or 'diagquadratic'.
%                                 Default: 'pseudolinear'
%       'Gamma'                 - Parameter for regularizing the
%                                 correlation matrix of predictors. For
%                                 linear discriminant, you can pass a
%                                 scalar greater or equal to 0 and less or
%                                 equal to 1. For quadratic discriminant,
%                                 you can pass either 0 or 1. Type 'help
%                                 fitcdiscr' for more info about the
%                                 'Gamma' parameter. Default: 0
%       'Delta'                 - Threshold on linear coefficients, a
%                                 non-negative scalar. For quadratic
%                                 discriminant, this parameter must be set
%                                 to 0. Default: 0
%       'FillCoeffs'            - If 'on', fills the struct holding
%                                 discriminant coefficients. If 'off', the
%                                 Coeffs property of discriminant DISCR
%                                 will be empty. Default: 'off'
%       'SaveMemory'            - If 'on', defers computing the full
%                                 covariance matrix until it is needed for
%                                 prediction. If 'off', computes and stores
%                                 the full covariance matrix in the
%                                 returned object. Set this parameter to
%                                 'on' if X has thousands of predictors.
%                                 Default: 'off'
%
%   See also ClassificationDiscriminant, fitensemble, fitcecoc, fitcdiscr.

%   Copyright 2013-2014 The MathWorks, Inc.
            
classreg.learning.FitTemplate.catchType(varargin{:});
temp = classreg.learning.FitTemplate.make('Discriminant','type','classification',varargin{:});
end
