function this = makecdiscr(mu,sigma,varargin)
%MAKECDISCR Make a discriminant from class means and covariance matrix.
%   DISCR=MAKECDISCR(MU,SIGMA) returns a compact discriminant for class
%   means MU and covariance matrix SIGMA. Use this method to construct a
%   compact discriminant from MU and SIGMA computed by tools other than
%   FITCDISCR. Use ClassificationDiscriminant/compact method to compact the
%   discriminant returned by FITCDISCR.
%
%   MU must be a K-by-P matrix of class means with one row per class and
%   one column per predictor. SIGMA must be either a P-by-P matrix or a
%   P-by-P-by-K array for P predictors and K classes. If SIGMA is P-by-P,
%   DISCR is a linear discriminant with pooled-in covariance matrix SIGMA.
%   If SIGMA is P-by-P-by-K, DISCR is a quadratic discriminant and
%   SIGMA(:,:,I) must return the covariance matrix for class I. The
%   pooled-in covariance matrix or the covariance matrix for each class
%   must be symmetric and positive semi-definite.
%
%   DISCR=MAKECDISCR(MU,SIGMA,'PARAM1',val1,'PARAM2',val2,...) specifies
%   optional parameter name/value pairs:
%       'BetweenSigma'     - Symmetric positive semi-definite matrix of
%                            size P-by-P with the between-class covariance
%                            values. If you pass 'BetweenSigma', its value
%                            is stored in the BetweenSigma property of
%                            DISCR. If you do not pass 'BetweenSigma', the
%                            BetweenSigma property is computed from the
%                            class means MU assuming equal weights for the
%                            classes. Default: []
%       'ClassNames'       - Array of class names with K elements.
%                            ClassNames can be a categorical predictor,
%                            character array, logical vector, numeric
%                            vector, or cell array of strings. If
%                            ClassNames is a character array, it must have
%                            one class label per row. Default: 1:K
%       'Cost'             - Square matrix, where COST(I,J) is the
%                            cost of classifying a point into class J if
%                            its true class is I. Alternatively, COST can
%                            be a structure S having two fields:
%                            S.ClassificationCosts containing the cost
%                            matrix C, and S.ClassNames containing the
%                            class names and defining the ordering of
%                            classes used for the rows and columns of the
%                            cost matrix. For S.ClassNames use the same
%                            data type as for 'ClassNames'. As an
%                            alternative, you can assign to the Cost
%                            property of discriminant DISCR. Default:
%                            COST(I,J)=1 if I~=J, and COST(I,J)=0 if I=J.
%       'FillCoeffs'       - If 'on', fills the struct holding discriminant
%                            coefficients. If 'off', the Coeffs property of
%                            discriminant DISCR will be empty. Default: 'on'
%       'PredictorNames'   - A cell array of names for the predictor
%                            variables with P elements, in the order in
%                            which they appear in MU. Default:
%                            {'x1','x2',...}
%       'Prior'            - Prior probabilities for each class.
%                               Specify as one of:
%                                   * String 'uniform' to set all class
%                                     probabilities equal
%                                   * A vector (one scalar value for each
%                                     class)
%                                   * A structure S with two fields:
%                                     S.ClassProbs containing a vector of
%                                     class probabilities, and S.ClassNames
%                                     containing the class names and
%                                     defining the ordering of classes used
%                                     for the elements of this vector.
%                               If you pass numeric values, MAKECDISCR
%                               normalizes them to add up to one. As an
%                               alternative, you can assign to the Prior
%                               property of discriminant DISCR. Default:
%                               'uniform'
%       'ResponseName'     - Name of the response variable, a string.
%                            Default: 'Y'
%
%   See also ClassificationDiscriminant, fitcdiscr,
%   ClassificationDiscriminant/compact.

%   Copyright 2013-2014 The MathWorks, Inc.

this = ClassificationDiscriminant.make(mu,sigma,varargin{:});
end
