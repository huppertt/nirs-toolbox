classdef kaiser < sigwin.parameterizewin
  %KAISER Construct a Kaiser object
  %
  %   SIGWIN.KAISER is not recommended. Use <a href="matlab:help kaiser">kaiser</a> instead.
  %
  %   H = SIGWIN.KAISER(N, B) constructs a Kaiser window object with length N
  %   and Beta B.  If N or B is not specified, they default to 64 and .5
  %   respectively.
  
  
  
  methods  % constructor block
    function hWIN = kaiser(n, param)
      
      hWIN.Name = 'Kaiser';
      createdynamicprops(hWIN, 'Beta', 'double','Beta');
      
      if nargin>0 && isnumeric(n),
        hWIN.Length = n;
      end
      
      if nargin>1,
        hWIN.Beta = param;
      else
        hWIN.Beta = 0.5;
      end
      
      
    end  % kaiser
    
  end  % constructor block
  
  methods  %% public methods
    function data=generate(hWIN)
      %GENERATE(hWIN) Generates the Kaiser window
      %
      %   sigwin.kaiser is not recommended.
      %   Use <a href="matlab:help kaiser">kaiser</a> instead.
      
      data = kaiser(hWIN.Length, hWIN.Beta);
    end
    
    function flag = isminordersupported(this) %#ok
      %ISMINORDERSUPPORTED Overloads the window base class method
      
      
      flag = 1;
    end
    
  end  %% public methods
  
end  % classdef

