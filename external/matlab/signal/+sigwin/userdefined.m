classdef (CaseInsensitiveProperties=true) userdefined < sigwin.window 
  %sigwin.userdefined class
  %   sigwin.userdefined extends sigwin.window.
  %
  %    sigwin.userdefined properties:
  %       Name - Property is of type 'string' (read only)
  %       MATLABExpression - Property is of type 'string'
  %
  %    sigwin.userdefined methods:
  %       generate - hWIN) Generates the userdefined window
  %       thisinfo - Information for this class.
  
  
  properties (AbortSet, SetObservable, GetObservable)
    %MATLABEXPRESSION Property is of type 'string'
    MATLABExpression = '';
  end
  
  properties (AbortSet, SetObservable, GetObservable, Hidden)
    %MATLAB_EXPRESSION Property is of type 'string' (hidden)
    MATLAB_expression = '';
  end
  
  
  methods  % constructor block
    function hWIN = userdefined(expression)
      %USERDEFINED Constructor of the userdefined class
      
      %   Author(s): V.Pellissier
      
      hWIN.Name = 'User Defined';
      if nargin > 0, hWIN.MATLABExpression = expression; end
      
      
    end  % userdefined
    
  end  % constructor block
  
  methods
    function set.MATLABExpression(obj,value)
      % DataType = 'string'
      validateattributes(value,{'char'}, {'row'},'','MATLABExpression')
      obj.MATLABExpression = value;
    end
    
    function value = get.MATLAB_expression(obj)
      value = getmatlab_expression(obj,obj.MATLAB_expression);
    end
    function set.MATLAB_expression(obj,value)
      % DataType = 'string'
      validateattributes(value,{'char'}, {'row'},'','MATLAB_expression')
      obj.MATLAB_expression = setmatlab_expression(obj,value);
    end
    
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the userdefined window
      
      %   Author(s): V.Pellissier
      %   Copyright 1988-2002 The MathWorks, Inc.
      
      try
        data = evalin('base', hWIN.MATLAB_expression);
      catch
        error(message('signal:sigwin:userdefined:generate:InvalidParam'));
      end
    end
    
    function [p, v] = thisinfo(h)
      %THISINFO Information for this class.
      
      % This should be a private method.
      
      %   Author(s): P. Costa
      %   Copyright 1988-2004 The MathWorks, Inc.
      
      p = {getString(message('signal:dfilt:info:MATLABExpression'))};
      v = {get(h, 'MATLAB_expression')};
    end
    
  end  %% public methods
  
end  % classdef

function me = setmatlab_expression(this, me)

set(this, 'MATLABExpression', me);
end  % setmatlab_expression


% -------------------------------------------------------------------------
function me = getmatlab_expression(this, me) %#ok

me = get(this, 'MATLABExpression');
end  % getmatlab_expression


% [EOF]
