function r = getblockinputprocessingrestrictions(~)
% GETBLOCKINPUTPROCESSINGRESTRICTIONS Get input processing restrictions for
% the block method.

%   Copyright 2012 The MathWorks, Inc.

% Does not support elements as channels
r = 'elementsaschannels';