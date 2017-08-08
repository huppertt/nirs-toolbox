function r = getrealizemdlraterestrictions(~,~)
% GETREALIZEMDLRATERESTRICTIONS Get rate restrictions for the realizemdl
% method.

%   Copyright 2012 The MathWorks, Inc.

% Single rate filters do not support any type of rate options
r = {'enforcesinglerate', 'allowmultirate'};