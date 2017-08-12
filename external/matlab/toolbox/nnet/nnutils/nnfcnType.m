classdef nnfcnType < nnfcnInfo
%NNTYPEFCNINFO Type function info.

% Copyright 2010 The MathWorks, Inc.

  properties (SetAccess = private)
    examples = {};
  end
  
  methods
    
    function x = nnfcnType(mname,name,version,examples)
      if nargin < 3, error(message('nnet:Args:NotEnough')); end
      if nargin < 4, examples = {}; end
      
      x = x@nnfcnInfo(['nntype.' mname],name,'nntype.type_fcn',version);
      x.examples = examples;
    end
    
    function disp(x)
      disp@nnfcnInfo(x)
    end
    
  end
  
end
