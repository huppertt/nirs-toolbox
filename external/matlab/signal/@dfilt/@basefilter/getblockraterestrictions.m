function r = getblockraterestrictions(~,inputProcessing)  %#ok<INUSD>
% GETBLOCKRATERESTRICTIONS Get rate restrictions for the block method.

%   Copyright 2012 The MathWorks, Inc.

% Single rate filters do not support any type of rate options
r = {'enforcesinglerate', 'allowmultirate'};
