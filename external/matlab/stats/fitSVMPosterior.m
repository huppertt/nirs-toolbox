function [obj,transformation] = fitSVMPosterior(obj,varargin)
%FITSVMPOSTERIOR Fit posterior probabilities for a Support Vector Machine model
%   OBJ=FITSVMPOSTERIOR(OBJ) finds the optimal transformation from scores
%   to posterior probabilities for an SVM model saved in object OBJ of type
%   ClassificationSVM or a cross-validated SVM model saved in object OBJ of
%   type ClassificationPartitionedModel. If you pass OBJ of type
%   ClassificationSVM, FITSVMPOSTERIOR computes SVM scores by 10-fold
%   cross-validation.
%
%   OBJ=FITSVMPOSTERIOR(OBJ,X,Y) finds the optimal transformation for an
%   SVM model saved in object OBJ of type CompactClassificationSVM using
%   predictor matrix X and class labels Y. Pass X as a floating-point
%   matrix with N rows. Pass Y as a character matrix with N rows or a
%   vector with N elements of one of the following types: single, double,
%   logical, categorical, and cell array of strings.
%
%   FITSVMPOSTERIOR returns an object of the same type as the input object
%   OBJ and sets the ScoreTransform property of this object to the optimal
%   transformation.
%
%   [OBJ,TRANS]=FITSVMPOSTERIOR(OBJ) or
%   [OBJ,TRANS]=FITSVMPOSTERIOR(OBJ,X,Y) also returns TRANS, a struct with
%   parameters of the optimal transformation from score S to posterior
%   probability P for the positive class in OBJ.ClassNames(2). TRANS has
%   the following fields:
%       Type                           - String, one of: 'sigmoid', 'step'
%                                        or 'constant'. If the two classes
%                                        overlap, FITSVMPOSTERIOR sets Type
%                                        to 'sigmoid'. If the two classes
%                                        are perfectly separated,
%                                        FITSVMPOSTERIOR sets Type to
%                                        'step'. If one of the two classes
%                                        has zero probability,
%                                        FITSVMPOSTERIOR sets Type to
%                                        'constant.
%          If Type is 'sigmoid', TRANS has additional fields:
%             Slope                    - Slope A of the sigmoid
%                                        transformation
%                                        P(S)=1/(1+exp(A*S+B))
%             Intercept                - Intercept B of the sigmoid
%                                        transformation
%                                        P(S)=1/(1+exp(A*S+B))
%          If Type is 'step', TRANS has additional fields:
%             PositiveClassProbability - Probability of the positive class
%                                        in the interval between LowerBound
%                                        and UpperBound
%             LowerBound               - Lower bound of the interval in
%                                        which the probability for the
%                                        positive class is set to
%                                        PositiveClassProbability. Below
%                                        this bound, the probability for
%                                        the positive class is zero.
%             UpperBound               - Upper bound of the interval in
%                                        which the probability for the
%                                        positive class is set to
%                                        PositiveClassProbability. Above
%                                        this bound, the probability for
%                                        the positive class is one.
%          If Type is 'constant', TRANS has additional fields:
%             PredictedClass           - Name of the predicted class, same
%                                        type as OBJ.ClassNames. The
%                                        posterior probability is one for
%                                        this class.
%
%   If OBJ is of type ClassificationSVM,
%   OBJ=FITSVMPOSTERIOR(OBJ,PARAM1',val1,'PARAM2',val2,...) specifies
%   optional parameter name/value pairs: 
%       'CVPartition'           - A partition created with CVPARTITION to
%                                 use for cross-validation. 
%       'Holdout'               - Holdout validation uses the specified
%                                 fraction of the data for test, and uses
%                                 the rest of the data for training.
%                                 Specify a numeric scalar between 0 and 1.
%       'KFold'                 - Number of folds to use for
%                                 cross-validation, a positive integer.
%                                 Default: 10
%       'Leaveout'              - Use leave-one-out cross-validation by
%                                 setting to 'on'.
%
%   See also ClassificationSVM,
%   classreg.learning.classif.CompactClassificationSVM,
%   classreg.learning.partition.ClassificationPartitionedModel, fitcsvm.

%   Copyright 2013-2014 The MathWorks, Inc.

if     isa(obj,'classreg.learning.classif.CompactClassificationSVM')
    
    obj = checkScoreTransform(obj);

    if isa(obj,'ClassificationSVM') && ...
            (isempty(varargin) || ischar(varargin{1}))
        cv = crossval(obj,varargin{:});
        Y = cv.Y;
        classnames = classreg.learning.internal.ClassLabel(cv.ClassNames);
        prior = cv.Prior;
        [~,s] = kfoldPredict(cv);        
        clear cv;
    else
        if numel(varargin)<2
            error(message('stats:fitSVMPosterior:PassXYtoCompactObject'));
        end
        
        X = varargin{1};
        Y = varargin{2};
        extraArgs = varargin(3:end);

        Y = classreg.learning.internal.ClassLabel(Y);
        
        if numel(Y)~=size(X,1)
            error(message('stats:fitSVMPosterior:XandYSizeMismatch'));
        end
        
        if numel(levels(Y))>2
            error(message('stats:fitSVMPosterior:TooManyLevelsInY'));
        end
        
        classnames = obj.ClassSummary.ClassNames;
        prior = obj.Prior;
        [~,s] = predict(obj,X,extraArgs{:});
    end
    
elseif isa(obj,'classreg.learning.partition.ClassificationPartitionedModel') ...
        && numel(obj.Trained)>0 ...
        && isa(obj.Trained{1},'classreg.learning.classif.CompactClassificationSVM')
    
    obj = checkScoreTransform(obj);
    
    if ~isempty(varargin)
        error(message('stats:fitSVMPosterior:DoNotPassOptionalParamsToCrossValidatedObject'));
    end
    
    Y = obj.Y;
    classnames = classreg.learning.internal.ClassLabel(obj.ClassNames);
    prior = obj.Prior;
    [~,s] = kfoldPredict(obj);
else
    error(message('stats:fitSVMPosterior:BadObjectType'));
end

if numel(classnames)~=2
    error(message('stats:fitSVMPosterior:ClassNamesMustHaveTwoElements'));
end

s = s(:,2);

Y = classreg.learning.internal.ClassLabel(Y);
y = grp2idx(Y,classnames) - 1;

if all(isnan(s))
    error(message('stats:fitSVMPosterior:AllScoresAreNaNs'));
end

% Return P=1 for negative class if positive class has zero probability
if sum(y==1)==0
    transformation.Type = 'constant';
    transformation.PredictedClass = labels(classnames(1));
    obj.ScoreTransform = eval(sprintf('@(S)constant(S,%i)',1));
    obj.ScoreType = 'probability';
    return;
end

% Return P=1 for positive class if negative class has zero probability
if sum(y==0)==0
    transformation.Type = 'constant';
    transformation.PredictedClass = labels(classnames(2));
    obj.ScoreTransform = eval(sprintf('@(S)constant(S,%i)',2));
    obj.ScoreType = 'probability';
    return;
end

smax0 = max(s(y==0));
smin1 = min(s(y==1));

if     smax0<=smin1 % perfect separation
    warning(message('stats:fitSVMPosterior:PerfectSeparation'));        
    transformation.Type = 'step';
    transformation.LowerBound = smax0;
    transformation.UpperBound = smin1;
    transformation.PositiveClassProbability = prior(2);
    f = eval(sprintf('@(S)step(S,%e,%e,%e)',smax0,smin1,prior(2)));
else    
    coeff = glmfit(s,y,'binomial','link','logit');
    a = -coeff(2);
    b = -coeff(1);
    transformation.Type = 'sigmoid';
    transformation.Slope = a;
    transformation.Intercept = b;
    f = eval(sprintf('@(S)sigmoid(S,%e,%e)',a,b));
end

obj.ScoreTransform = f;
obj.ScoreType = 'probability';

end

function out = constant(in,cls) %#ok<DEFNU>
out = zeros(size(in));
out(:,cls) = 1;
end

function out = sigmoid(in,a,b) %#ok<DEFNU>
out = zeros(size(in));
out(:,2) = 1./(1+exp(a*in(:,2)+b));
out(:,1) = 1 - out(:,2);
end

function out = step(in,lo,hi,p) %#ok<DEFNU>
out = zeros(size(in));
s = in(:,2);
out(s>hi,2) = 1;
out(s<lo,1) = 1;
between = s>=lo & s<=hi;
out(between,2) = p; 
out(between,1) = 1-p; 
end

function obj = checkScoreTransform(obj)
if ~strcmp(obj.ScoreTransform,'none')
    warning(message('stats:fitSVMPosterior:ResetScoreTransform'));
    obj.ScoreTransform = 'none';
end
end
