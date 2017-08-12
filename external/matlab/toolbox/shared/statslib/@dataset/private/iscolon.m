function tf = iscolon(indices)
%ISCOLON Check if a set of indices is ':'.

%   Copyright 2007 The MathWorks, Inc. 


% Check ischar first.  isequal(58,':') alone is true,
% and strcmp({':'},':') is also true
tf = ischar(indices) && strcmp(indices,':');
