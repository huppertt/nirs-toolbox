classdef CompactClassifByBinaryRegr < classreg.learning.classif.ClassificationModel
%CompactClassifByBinaryRegr Binary classification by regression.

%   Copyright 2010-2014 The MathWorks, Inc.

   
    properties(GetAccess=public,SetAccess=protected)
        CompactRegressionLearner = []; % Compact regression model
    end
    
    methods(Access=protected)        
        function this = CompactClassifByBinaryRegr(...
                dataSummary,classSummary,scoreTransform,crl)
            this = this@classreg.learning.classif.ClassificationModel(...
                dataSummary,classSummary,scoreTransform,[]);
            this.CompactRegressionLearner = crl;
        end
        
        function s = score(this,X,varargin)
            s = predict(this.CompactRegressionLearner,X,varargin{:});
            if numel(this.ClassSummary.ClassNames)>1
                s = [s(:) -s(:)];
            end
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.ClassificationModel(this,s);
            s.CompactRegressionLearner = this.CompactRegressionLearner;
        end
    end

end
