function funloss = lossCheck(funloss,type)  

%   Copyright 2010-2014 The MathWorks, Inc.


if     ischar(funloss)
    if     strcmp(type,'classification')
        allowed = {'binodeviance' 'classifedge' 'classiferror' ...
            'exponential' 'mincost' 'hinge' 'quadratic'};
    elseif strcmp(type,'regression')
        allowed = {'mse'};
    else
        allowed = {};
    end
    idx = find(strncmpi(funloss,allowed,length(funloss)));
    if isempty(idx) || ~isscalar(idx)
        error(message('stats:classreg:learning:internal:lossCheck:BadFunlossString'));
    end
    funloss = str2func(['classreg.learning.loss.' allowed{idx}]);
elseif ~isa(funloss,'function_handle')
    error(message('stats:classreg:learning:internal:lossCheck:BadFunlossType'));
end
end
