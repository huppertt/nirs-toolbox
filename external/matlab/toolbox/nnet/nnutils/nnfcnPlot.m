classdef nnfcnPlot < nnfcnInfo
%NNPLOTFCNINFO Plot function info.

% Copyright 2010 The MathWorks, Inc.

  properties (SetAccess = private)
  end
  
  methods
    
    function x = nnfcnPlot(name,title,version,param)
      if nargin < 3, error(message('nnet:Args:NotEnough')); end
      
      x = x@nnfcnInfo(name,title,'nntype.plot_fcn',version);
      x.setupParameters(param);
    end
    
    function disp(x)
      disp@nnfcnInfo(x)
      fprintf('\n')
      disp(' <a href="matlab:doc nnPlotFunctionInfo">nnPlotFunctionInfo</a>')
      fprintf('\n')
      %=======================:
      disp('           (none) ');
    end
    
  end
  
end

