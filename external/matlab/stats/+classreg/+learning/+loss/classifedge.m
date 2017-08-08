function e = classifedge(C,Sfit,W,cost)

%   Copyright 2010 The MathWorks, Inc.



m = classreg.learning.loss.classifmargin(C,Sfit);

% Average the margins
notNaN = ~isnan(m);
e = sum(W(notNaN).*m(notNaN)) / sum(W(notNaN));
end
