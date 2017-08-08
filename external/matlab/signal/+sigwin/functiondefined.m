classdef (CaseInsensitiveProperties=true) functiondefined < sigwin.variablelength 
  %sigwin.functiondefined class
  %   sigwin.functiondefined extends sigwin.variablelength.
  %
  %    sigwin.functiondefined properties:
  %       Name - Property is of type 'string' (read only)
  %       Length - Property is of type 'spt_uint32 user-defined'
  %       MATLABExpression - Property is of type 'string'
  %       Parameters - Property is of type 'mxArray'
  %
  %    sigwin.functiondefined methods:
  %       generate - hWIN) Generates the functiondefined window
  
  
  properties (AbortSet, SetObservable, GetObservable)
    %MATLABEXPRESSION Property is of type 'string'
    MATLABExpression = '';
    %PARAMETERS Property is of type 'mxArray'
    Parameters = [];
  end
  
  properties (AbortSet, SetObservable, GetObservable, Hidden)
    %MATLAB_EXPRESSION Property is of type 'string' (hidden)
    MATLAB_expression = '';
  end
  
  
  methods  % constructor block
    function hWIN = functiondefined(fcnname, n, params)
      %FDEFWIN Constructor of the functiondefined class
      
      %   Author(s): V.Pellissier
      
      % hWIN = sigwin.functiondefined;
      hWIN.Name = 'User Defined';
      
      if nargin>0,
        hWIN.MATLABExpression = fcnname;
      end
      
      if nargin>1 && isnumeric(n),
        hWIN.Length = n;
      end
      
      if nargin>2,
        hWIN.Parameters = params;
      end
      
      
    end  % functiondefined
    
    
    function set.MATLABExpression(obj,value)
      % DataType = 'string'
      validateattributes(value,{'char'}, {'row'},'','MATLABExpression')
      obj.MATLABExpression = value;
    end
    
    function set.Parameters(obj,value)
      obj.Parameters = value;
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
      %GENERATE(hWIN) Generates the functiondefined window
      
      if isempty(hWIN.Parameters),
        try
          data = feval(hWIN.MATLAB_expression, hWIN.Length);
        catch ME
          throw(ME);
        end
      else
        params = hWIN.Parameters;
        if ~iscell(params),
          params = {params};
        end
        try
          data = feval(hWIN.MATLAB_expression, hWIN.Length, params{:});
        catch ME
          throw(ME);
        end
      end
    end
    
  end  %% public methods
  
end  % classdef

function me = setmatlab_expression(this, me)

set(this, 'MATLABExpression', me);
me = '';
end  % setmatlab_expression


% -------------------------------------------------------------------------
function me = getmatlab_expression(this, me) %#ok

me = get(this, 'MATLABExpression');
end  % getmatlab_expression


% [EOF]
