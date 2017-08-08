classdef (CaseInsensitiveProperties=true) chebwin < sigwin.parameterizewin
  %CHEBWIN Construct a Chebyshev object
  %
  %   SIGWIN.CHEBWIN is not recommended.  Use <a href="matlab:help chebwin">chebwin</a> instead.
  %
  %   H = SIGWIN.CHEBWIN(N, S) constructs a Chebyshev window object with
  %   length N and sidelobe attenuation S.  If N or S is not specified, they
  %   default to 64 and 100 respectively.
  
  
  
  properties (Transient, SetObservable, GetObservable, Hidden)
    %SIDELOBE_ATTEN Property is of type 'double' (hidden)
    Sidelobe_atten
  end
  
  
  methods  % constructor block
    function hWIN = chebwin(n, atten)
      
      
      % hWIN = sigwin.chebwin;
      hWIN.Name = 'Chebyshev';
      createdynamicprops(hWIN, 'SidelobeAtten', 'double', 'Sidelobe Attenuation');
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      if nargin>1,
        hWIN.SidelobeAtten = atten;
      else
        hWIN.SidelobeAtten = 100;
      end
      
      
    end  % chebwin
    
    
    function value = get.Sidelobe_atten(obj)
      value = getsidelobe_atten(obj,obj.Sidelobe_atten);
    end
    function set.Sidelobe_atten(obj,value)
      % DataType = 'double'
      validateattributes(value,{'double'}, {'scalar'},'','Sidelobe_atten')
      obj.Sidelobe_atten = setsidelobe_atten(obj,value);
    end
    
    
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Chebyshev window
      %
      %   sigwin.chebwin is not recommended.
      %   Use <a href="matlab:help chebwin">chebwin</a> instead.
      
      
      data = chebwin(hWIN.Length, hWIN.SidelobeAtten);
      
    end
    
  end  %% public methods
  
end  % classdef

function sa = setsidelobe_atten(this, sa) %#ok<INUSD>

error(message('signal:sigwin:chebwin:schema:DeprecatedProperty'));
end  % setsidelobe_atten


% -------------------------------------------------------------------------
function sa = getsidelobe_atten(this, sa) %#ok

error(message('signal:sigwin:chebwin:schema:DeprecatedProperty'));
end  % getsidelobe_atten



% [EOF]
