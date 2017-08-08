function [outputs,flag] = finddesignfiltflag(~,inputs)
%FINDDESIGNFILTFLAG   Find FromDesignfilt flag

%   Copyright 2013 The MathWorks, Inc.

flag = false;
for idx = 1:numel(inputs)
  currentInput = inputs{idx};
  if ischar(currentInput) && strcmpi(currentInput,'FromDesignfilt')
    inputs(idx) = [];
    flag = true;
  end
end
   
outputs = inputs;

