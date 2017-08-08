function temp = templateECOC(varargin)
%TEMPLATEECOC Create a template for ECOC learning.
%   T=TEMPLATEECOC() returns a template suitable for multiclass learning.
%
%   T=TEMPLATEECOC('PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%       'Coding'        - Either string or matrix of size K-by-L for K classes
%                         and L binary learners.
%                           * If you pass 'Coding' as a string, pass it as one
%                             of:
%                               'onevsone' or 'allpairs'
%                               'onevsall'
%                               'binarycomplete'
%                               'ternarycomplete'
%                               'ordinal'
%                               'sparserandom'
%                               'denserandom'
%                             See help for DESIGNECOC for details.
%                           * If you pass 'Coding' as a matrix, compose this
%                             matrix using the following rules:
%                               - Every element of this matrix must be one of:
%                                 -1, 0, and +1.
%                               - Every column must contain at least one -1 and
%                                 at least one +1.
%                               - There are no equal columns or columns equal
%                                 after a sign flip.
%                               - There are no equal rows.
%                             Default: 'onevsone'
%       'FitPosterior'  - True or false. If true, classification scores returned
%                         by the binary learners are transformed to posterior
%                         probabilities and the PREDICT method of the ECOC model
%                         OBJ can compute ECOC posterior probabilities. Default:
%                         false
%       'Learners'      - A cell array of binary learner templates or a single
%                         binary learner template. You must construct every
%                         template by calling one of: templateDiscriminant,
%                         templateEnsemble, templateKNN, templateSVM, and
%                         templateTree. If you pass 'Learners' as a single
%                         template object, FITCECOC will construct all binary
%                         learners using this template. If you pass a cell array
%                         of objects, its length must match the number of
%                         columns in the coding matrix. FITCECOC then uses the
%                         l-th template to construct a binary learner for the
%                         l-th column of the coding matrix. If you pass one
%                         binary learner with default parameters, you can pass
%                         'Learners' as a string with the name of the binary
%                         learner, for example, 'SVM'. Default: 'SVM'
%
%   See also ClassificationECOC, fitcecoc, designecoc, templateDiscriminant,
%   templateEnsemble, templateKNN, templateSVM, templateTree.

%   Copyright 2014 The MathWorks, Inc.

temp = classreg.learning.FitTemplate.make('ECOC',varargin{:});
end
